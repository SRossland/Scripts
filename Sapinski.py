# This is a program to do the sapinkski triangle.

# Rules: 3 points must be made in such a configuration that the terminus of 
# each point represents the point of an equilatiral triangle.

# Pick a random point that is the halfway point between two points. This new point
# will be one vertex of a new triangle that you will repeat this last step.


import turtle
import random

screen = turtle.Screen()
screen.title('Triangle Chaos Game with Python Turtle')
screen.setup(1000,1000)
screen.tracer(0,0)
turtle.hideturtle()
turtle.speed(0)
turtle.up()

A = (0,350)
B = (-300,-200)
C = (300,-200)
V = (A,B,C) # list of three vertices
for v in V:
    turtle.goto(v)
    turtle.dot('red')

n = 10000 # number of points to draw
p = (random.uniform(-200,200),random.uniform(-200,200)) # random starting point
t = turtle.Turtle()
t.up()
t.hideturtle()
for i in range(n):
    t.goto(p)
    t.dot(2,'blue')
    r = random.randrange(len(V)) # pick a random vertex
    p = ((V[r][0]+p[0])/2,(V[r][1]+p[1])/2) # go to mid point between the random vertex and point   
    if i % 1000 == 0: # update for every 1000 moves, this part is for performance reason only
        t = turtle.Turtle() # use new turutle
        t.up()
        t.hideturtle()
        screen.update()






