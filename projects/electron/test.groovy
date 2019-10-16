x = [1, 2, 5, 6]
y1 = [false, true, true, true]
y2 = [true, false, false, false]

def everythingElsePassed(results, i) {

     for (j=0; j<results.size(); j++){
     	 if (i != j && results[j] == false){
	    return false 
	 }
     }
     

    return true
}

def everythingElseFailed(results, i){
    def reversed = results.collect{ !it }
    return everythingElsePassed(reversed,i)
}

println(everythingElsePassed(y1, 0))
println(everythingElsePassed(y1, 1))
println(everythingElsePassed(y1, 2))
println(everythingElsePassed(y1, 3))

println(everythingElseFailed(y2,0))
println(everythingElseFailed(y2,1))
println(everythingElseFailed(y2,2))
println(everythingElseFailed(y2,3))
