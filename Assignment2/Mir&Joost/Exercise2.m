%% SC4210 NDC: Assignment 1, Exercise 1
clear all; close all; clc;
W = [1 1 0 1 0 0
     1 0 1 0 1 0
     0 1 0 0 0 1
     1 0 0 0 1 0
     0 1 0 1 0 1
     0 0 1 0 1 0];
w = sum(W,2);
D = diag(w);
P = inv(D)*W;
L = D-W;


 