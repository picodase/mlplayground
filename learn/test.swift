print("hello world, but in Swift!")

let a = 67
var b = 9

let weight: Float = 5.0
print(weight)

var shopping_list = ["carrots", "beanpaste"]

var diction = [
    1:"hello",
    2:"world"
]

print(diction)

let listobisto=["hello", "yellow"]

for i in listobisto{
    print(i)
}

var listobistu=[0, 2]

for i in listobistu{
    print(i*2)
}

for i in 0...21{
    print(i)
}

// commentationationation

/*
multi-line comments, ohhhyeah!
*/

func yelloooooo(jaundice:Int) -> Bool {
    if jaundice==1{
        return true
    }
    else{
        return false
    }
}

yelloooooo(jaundice:1)

if (yelloooooo(jaundice:0)==false) {
    print("yes, they are equal!")
}
else{
    print("no, unequal")
}

import Python

// Load numpy from Python
let np = Python.import("numpy")

// Create an array of zeros
var zeros = np.ones([2, 3])
print(zeros)