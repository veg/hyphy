function foo (x) {
    Call ("bar", x);
}

function bar (y) {
    return "Beavis";
}

fprintf (stdout, Call ("foo", 5));

fprintf (stdout, "\nIn between\n");

fprintf (stdout, Call ("foo", 5));

/*

ONE

     [1] 1 -- create _FString
     [2] 2 -- return statement from Bar 
     [3] 3 -- bar.chain result 
     [4] 4 -- Call (bar) local Formula result
     [5] 5 -- +reference inside Call (bar)
     [4] 6 -- Deleting local Formula result (for Call/bar)
     [5] 7 -- Computing the Call (bar, x) inside foo
     [4] 8 -- Decompiling Call ("bar", x)
     [3] 9 -- clear stack result of the above function
     [2] 10-- decompiling return "Beavis"
     [1] 11-- removing return "Beavis" _Stack 
     [0] 12-- Deleting _FString
TWO     

     [1+] 1 -- create _FString
     [2+] 2 -- _Stack for the return statement from "bar" 
     [3+x] 3 -- bar.chain result 
     [4+x] 4 -- Call ("bar", x) local Formula result _Stack
     [5+] 5 -- +reference inside .Call
     [4-] 6 -- Deleting local Formula result (for Call/bar) / releases Step 4
     [5+x] 7 -- Computing  Call (bar, x) inside foo; returning via cachedResult
     [4-] 8 -- Clear stack of Call (bar, x) inside foo
     [3-] 9 -- Clear "bar" chain result // match 3
     [2-]10 -- Clear the _Stack of the bar return command // match 4 
     
     at this point the references are being held by the original _FString 
     AND the return statement from Bar?
     
     
     [3]11 -- Compute the "bar" return result
     [4]12 -- "bar" return + reference
     [5]13 -- "bar" Call "the_call" compute
     [6]14 -- "bar" Call result + ref
     [5]15 -- "bar" Call clear the_call
     [6]16 -- Computing the Call (bar, x) inside foo // compiled    
     [5]17 -- Decompiling Call ("bar", x)
     [4]18 -- Decompiling Call ("bar", x) / stack clear
     [3]17 -- Decompiling return from "bar"
     [2]18 -- Decompiling return from "bar" / stack clear
     [1]19 -- Deleting "result" from "bar"
    
*/