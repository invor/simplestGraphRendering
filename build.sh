#!/bin/bash

cd src

g++ -W -Wall --std=c++11 simplestGraphRendering.cpp -o simplestGraphRendering -lGL -lglut -lm -lX11 -lglfw -lGLEW

./simplestGraphRendering -gf ../resources/test.gl 
