/* The Computer Language Benchmarks Game
   http://benchmarksgame.alioth.debian.org/

   contributed by Brad Chamberlain

   derived from the GNU C version by Аноним Легионов and Jeremy Zerfas
     as well as from previous Chapel versions by Casey Battaglino,
     Kyle Brady, and Preston Sahabu.
*/

// based on a cobegin-based fasta-blc

use List;

config const n = 1000,           // the length of the generated strings
             lineLength = 60,    // the number of columns in the output
             blockSize = 1024;   // the parallelization granularity

config const frames = 5;         // number of pipeline frames to store
var allocated_frames = 0;

const chunkSize = lineLength*blockSize;

config param debugFasta = false;

//
// Nucleotide definitions
//
enum nucleotide {
  A = ascii("A"), C = ascii("C"), G = ascii("G"), T = ascii("T"),
  a = ascii("a"), c = ascii("c"), g = ascii("g"), t = ascii("t"),
  B = ascii("B"), D = ascii("D"), H = ascii("H"), K = ascii("K"),
  M = ascii("M"), N = ascii("N"), R = ascii("R"), S = ascii("S"),
  V = ascii("V"), W = ascii("W"), Y = ascii("Y")
}
use nucleotide;

//
// Sequence to be repeated
//
const ALU: [0..286] int(8) = [
  G, G, C, C, G, G, G, C, G, C, G, G, T, G, G, C, T, C, A, C,
  G, C, C, T, G, T, A, A, T, C, C, C, A, G, C, A, C, T, T, T,
  G, G, G, A, G, G, C, C, G, A, G, G, C, G, G, G, C, G, G, A,
  T, C, A, C, C, T, G, A, G, G, T, C, A, G, G, A, G, T, T, C,
  G, A, G, A, C, C, A, G, C, C, T, G, G, C, C, A, A, C, A, T,
  G, G, T, G, A, A, A, C, C, C, C, G, T, C, T, C, T, A, C, T,
  A, A, A, A, A, T, A, C, A, A, A, A, A, T, T, A, G, C, C, G,
  G, G, C, G, T, G, G, T, G, G, C, G, C, G, C, G, C, C, T, G,
  T, A, A, T, C, C, C, A, G, C, T, A, C, T, C, G, G, G, A, G,
  G, C, T, G, A, G, G, C, A, G, G, A, G, A, A, T, C, G, C, T,
  T, G, A, A, C, C, C, G, G, G, A, G, G, C, G, G, A, G, G, T,
  T, G, C, A, G, T, G, A, G, C, C, G, A, G, A, T, C, G, C, G,
  C, C, A, C, T, G, C, A, C, T, C, C, A, G, C, C, T, G, G, G,
  C, G, A, C, A, G, A, G, C, G, A, G, A, C, T, C, C, G, T, C,
  T, C, A, A, A, A, A
];

//
// Index aliases for use with (nucleotide, probability) tuples
//
param nucl = 1,
      prob = 2;

param IM = 139968;

//
// Probability tables for sequences to be randomly generated
//
const IUB = [(a, 0.27), (c, 0.12), (g, 0.12),
                                     (t, 0.27), (B, 0.02), (D, 0.02),
                                     (H, 0.02), (K, 0.02), (M, 0.02),
                                     (N, 0.02), (R, 0.02), (S, 0.02),
                                     (V, 0.02), (W, 0.02), (Y, 0.02)];

const HomoSapiens = [(a, 0.3029549426680),
                                            (c, 0.1979883004921),
                                            (g, 0.1975473066391),
                                            (t, 0.3015094502008)];

//
// Redefine stdout to use lock-free binary I/O and capture a newline
//
const stdout = openfd(1).writer(kind=iokind.native, locking=false);
param newline = ascii("\n"): int(8);

class WorkChunk {
  var frame: int;
  var startIdx: int;
  var length: int;
  var rands: [0..chunkSize] int/*(32)*/;
  var line_buff: [0..(lineLength+1)*blockSize-1] int(8);
}

record WorkList {
  var llist: list(WorkChunk);
  var lock$: sync bool;
  var wait$: sync bool;
  var done: bool;
  var name: string;
  var ordered = false;
  var lastFrame = -1; // last frame returned
  var framesRemaining: int; // frames left before this list is done


  // XXX List.list doesn't support an "insert_after()", can only add
  // at the beginning or end.  So we can't maintain llist as an
  // ordered list.  Instead, append each new wc (unless its frame is <
  // the llist first item's frame, in which case prepend it.)
  // Then do more work at get() time.
  proc append(wc: WorkChunk) {
    // The WorkList lock is held, put this wc on llist.
    proc _append(wc: WorkChunk) {
      if ordered && llist.length > 0 {
	var atStart: bool;
	for first in llist {
	  if wc.frame < first.frame {
	    atStart = true;
	  } else {
	    atStart = false;
	  }
	  break;
	}
	if (atStart) {
	  llist.prepend(wc);
	} else {
	  llist.append(wc);
	}
      } else {
	llist.append(wc);
      }
    }
    var wasEmpty: bool;
    lock$ = true;
    debug(name, " append: got lock");
    //    debugAssert(!done, name, " append: !done"); doesn't apply to freeList now
    var needsWake = !wait$.isFull; // could be due to llist empty, or out-of-order
    _append(wc);
    if (needsWake) {
      debug(name, " append: clearing wait -> full");
      wait$ = false;
    }
    debug(name, " append: done ", wc.frame, " length ", llist.length);
    lock$;
  }
  
  proc get(blocking = true) {
    proc _get() {
      var wc: WorkChunk;
      if ordered {
	var i: WorkChunk;
	for i in llist {
	  if i.frame == lastFrame + 1 {
	    wc = i;
	    break;
	  }
	}
	if wc != nil {
	  llist.remove(wc);
	}
      } else {
	wc = llist.pop_front();
      }
      if llist.length == 0 || (wait$.isFull && wc == nil) {
	debug(name, " emptying wait$");
	debugAssert(wait$.isFull, name, " wait$ not full");
	wait$;
      }
      if wc != nil {
	lastFrame = wc.frame;
	framesRemaining -= 1;
	if (framesRemaining == 0) {
	  _setDone();
	}
	debug(name, " get returning ", wc.frame, ", remaining: ", framesRemaining);
      }
      return wc;
    }

    var wc: WorkChunk;
    while (true) {
      debug(name, " get loop");
      lock$ = true;
      debug(name, " get locked, length ", llist.length);
      if llist.length > 0 {
	wc = _get();
      }
      lock$;
      if wc {
	debug(name, " get got ", wc.frame);
	return wc;
      }
      if done {
	debug(name, " get got done");
	return nil;
      }
      if !blocking {
	debug(name, " get not blocking");
	return nil;
      }
      debug(name, " get waiting");
      wait$.readFF();
    }
    assert(false);
    return nil;
  }

  proc _setDone() {
    debug(name, " _setDone");
    done = true;
    debugAssert(wait$.isFull == (llist.length > 0), name, " isFull ", wait$.isFull, " and length ", llist.length);
    if llist.length == 0 {
      wait$ = false;
    }
  }

  proc reinit() {
    lock$ = true;
    debug(name, " reinit length", llist.length);
    done = false;
    framesRemaining = (n + chunkSize - 1) / chunkSize;
    if llist.length > 0 { // Should only be freeList
      debugAssert(wait$.isFull);
    } else {
      if (wait$.isFull) {
	wait$;
      }
    }
    if ordered {
      lastFrame = -1;
    }
    lock$;
  }
}

var freeList, randList, filledList: WorkList;
filledList.ordered = true;

proc main() {
  allocWorkChunks();
  
  repeatMake(">ONE Homo sapiens alu\n", ALU, 2*n);
  randomMake(">TWO IUB ambiguity codes\n", IUB, 3*n);
  randomMake(">THREE Homo sapiens frequency\n", HomoSapiens, 5*n);

  debug("needed to allocate ", allocated_frames, " extra frames");
}

proc allocWorkChunks() {
  if debugFasta {
    // Line up the names printed in debug messages.
    freeList.name = "freeList";
    randList.name = "                         randList";
    filledList.name = "                                                filledList";
  } else {
    freeList.name = "freeList";
    randList.name = "randList";
    filledList.name = "filledList";
  }
  
  for i in 1..frames {
    var wc = new WorkChunk();
    freeList.append(wc);
  }
}

proc reinitLists() {
  debug("reiniting lists");
  freeList.reinit();
  randList.reinit();
  filledList.reinit();
  debug("lists reinited");
}


//
// Repeat string 'str' for 'n' characters
//
proc repeatMake(desc, str, n) {
  stdout.write(desc);

  const r = str.size,
        s = [i in 0..(r+lineLength)] str[i % r];

  for i in 0..n by lineLength {
    const lo = i % r + 1,
          len = min(lineLength, n-i);
    stdout.write(s[lo..#len], newline);
  }
}

//
// Output a random sequence of length 'n' using distribution 'a'
//
proc randomMake(desc, nuclInfo, n) {
  const numNucls = nuclInfo.size;

  stdout.write(desc);
  stdout.flush();
  
  var cumul_p: [1..numNucls] int;
  //
  // Sum the probabilities of the nucleotide info
  //
  var p = 0.0;
  for i in 1..numNucls {
    p += nuclInfo[i](prob);
    cumul_p[i] = 1 + (p*IM):int;
  }


  if false {
  if (here.maxTaskPar < 3) then
    writef("Warning: This code uses busy waiting and 3 tasks, so may deadlock"+
           " since here.maxTaskPar = %i", here.maxTaskPar);
  }

  debugAssert(freeList.llist.length == frames + allocated_frames);
  reinitLists();

  debug("finished reinit");

  cobegin {
    computeRands();
    computeLines();
    computeLines();
    writeLines();
  }

  proc computeRands() {
    var frame = 0;
    debug("computeRands: 1..", n, " by ", chunkSize);
    for i in 1..n by chunkSize {
      const bytes = min(chunkSize, n-i+1);

      debug("computeRands: i ", i, " bytes ", bytes);
      var wc = freeList.get(blocking = false);
      if (wc == nil) {
	debug("allocating new wc");
	wc = new WorkChunk();
	allocated_frames += 1;
      }
      debugAssert(wc != nil);
      debug("got a wc for rand #", frame);
      
      wc.frame = frame;
      wc.startIdx = i;
      wc.length = bytes;

      getRands(bytes, wc.rands);

      randList.append(wc);

      frame += 1;
    }
  }

  proc computeLines() {
    var frame = 0;
    while true {
      var wc = randList.get();
      if (wc == nil) {
	break;
      }

      var bytes = wc.length;
      var col = 0;
      var off = 0;
      ref myRands = wc.rands,
	  myBuff = wc.line_buff;
      for i in 0..#bytes {
        const r = myRands[i];
        var ncnt = 1;
        for j in 1..numNucls do
          if r >= cumul_p[j] then
            ncnt += 1;

        myBuff[off] = nuclInfo[ncnt](nucl);

        off += 1;
        col += 1;
        if (col == lineLength) {
          col = 0;
          myBuff[off] = newline;
          off += 1;
        }
      }
      if (col != 0) {
        myBuff[off] = newline;
        off += 1;
      }

      wc.length = off;
      filledList.append(wc);
    }
}

  proc writeLines() {
    var frame = -1;
    while true {
      var wc = filledList.get();
      if (wc == nil) {
	return;
      }

      if (wc.frame != frame + 1) {
	assert(false, "Oops, got output frame ", wc.frame, " after frame ", frame);
      }

      var off = wc.length;
      stdout.write(wc.line_buff[0..#off]);

      frame = wc.frame;

      freeList.append(wc);
    }
  }
}


//
// Deterministic random number generator
//
var lastRand = 42/*:int(32)*/;

proc getRands(n, arr) {
  param IA = 3877,
        IC = 29573;

  //  writef("tid %i got turn\n", tid);
  for i in 0..#n {
    lastRand = (lastRand * IA + IC) % IM;
    arr[i] = lastRand;
  }
}



inline proc debug(x ... ?k) {
  if debugFasta {
    stderr.writeln((... x));
    stderr.flush();
  }
}

inline proc debugAssert(x ... ?k) {
  if debugFasta {
    assert((...x));
  }
}
