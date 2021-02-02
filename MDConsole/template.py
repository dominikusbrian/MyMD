#!/usr/bin/env python3


# Python templates for development purpose
class MyClass:
    def __init__(self,parm1, parm2):
        self.parm1 = parm1
        self.parm2 = parm2
    
    def myfunction(self):
        print("test,test, this is " + self.parm1 )

#Test Deploy zone
p1 = MyClass("parm1", 33)

print(p1.parm1)
print(p1.parm2)

p1.myfunction()
