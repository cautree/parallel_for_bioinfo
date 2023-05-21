println("Hello, World!")
my_var = 3
println(my_var)

chr = "1"
println("processing chromosome $chr")

kmers = [11,21,27,31]
println(kmers[0])


kmers = [11,21,27,31]
//Lists can also be indexed with negative indexes
println(kmers[3])
println(kmers[-1])


kmers = [11,21,27,31]
// The first three elements using a range.
println(kmers[0..2])


kmers = [11,21,27,31]
println("The first three elements in the Lists are. $kmers[0..2]")
println("The first three elements in the Lists are. ${kmers[0..2]}")


mylist = [0,1,2]
println(mylist.size())
//inside a string need we need to use the ${} syntax
println("list size is:  ${mylist.size()}")

mylist = [0,1,2]
println mylist.get(1)


mylist = [1,2,3]
println mylist
println mylist + [1]
println mylist - [1]
println mylist * 2
println mylist.reverse()
println mylist.collect{ it+3 }
println mylist.unique().size()
println mylist.count(1)
println mylist.min()
println mylist.max()
println mylist.sum()
println mylist.sort()
println mylist.find{it%2 == 0}
println mylist.findAll{it%2 == 0}



list = [1,2,3,4,5,6,7,8,9,10]
//or
list = 1..10
println("${list[4]}")
//or
println("${list.get(4)}")


roi = [ chromosome : "chr17", start: 7640755, end: 7718054, genes: ['ATP1B2','TP53','WRAP53']]

//Use of the square brackets.
println(roi['chromosome'])

//Use a dot notation            
println(roi.start)

//Use of get method                      
println(roi.get('genes'))

//Use of the square brackets
roi['chromosome'] = '17'

//Use a dot notation        
roi.chromosome = 'chr17'  

//Use of put method              
roi.put('genome', 'hg38')  


square = { it * it }

square = { it * it }
x = [ 1, 2, 3, 4 ]
y = x.collect(square)
println y


x = [ 1, 2, 3, 4 ]
y = x.collect({ it * it })
println("x is $x")
println("y is $y")

x=[1,2,3,4,5,6]
chrx = { x -> "chr" + x}
y = x.collect(chrx)
println(y)

prefix = { "chr${it}"}
y = x.collect(prefix)
print(y)


tp53 = [chromosome: "chr17",start:7661779 ,end:7687538, genome:'GRCh38', gene: "TP53"]
//perform subtraction of end and start coordinates
region_length = {start,end -> end-start }
tp53.length = region_length(tp53.start,tp53.end)
println(tp53)



//closure with two parameters
printMap = { a, b -> println "$a with value $b" }

//map object
my_map = [ chromosome : "chr17", start : 1, end : 83257441 ]

//each iterates through each element
my_map.each(printMap)



x = 12
if( x > 10 )
    println "$x is greater than 10"
    
    
list = [1,2,3]
if( list )
    println list
else
    println 'The list is empty'
    
println list ? list : 'The list is empty'


for (int i = 0; i <3; i++) {
   println("Hello World $i")
}

list = ['a','b','c']

for( String elem : list ) {
  println elem
}


int fib(int n) {
    return n < 2 ? 1 : fib(n-1) + fib(n-2)
}

println (fib(10)) // prints 89



def fact( n ) {
  n > 1 ? n * fact(n-1) : 1
}

println (fact(5)) // prints 120

