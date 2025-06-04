clc
clear all
addpath('../ROUTINES/QUADRATURE/')

%% Gauss-Legendre Quadrature (accuracy up to 2n-1)
% Represents:
% \int_{-1}^1 f(x) dx = \sum_i wi f(xi)
Nq = 5;
[xi, wi] = LGWT(Nq);

[pl, intg] = PLEGE(0, xi);

[wi'*(pl.^2) intg]

%% Gauss-Lobatto Quadrature (end points included; accuracy is up to 2n-3)
% Represents:
% \int_{-1}^1 f(x) dx = \sum_i wi f(xi)
Nq = 5;
[xi, wi] = LGLWT(Nq);

[pl, intg] = PLEGE(2, xi);

[wi'*(pl.^2) intg]

%% Probabilists' Gauss-Hermite Quadrature (accurate up to 2n-1)
% Represents:
% \int_{\inf}^{\inf} f(x) \frac{1}{\sqrt{2\pi}}exp(-x^2/2) dx = \sum_i wi f(xi)

Nq = 5;
[xi, wi] = GPHWT(Nq);

[pl, intg] = PHERM(0, xi);

[wi'*(pl.^2) intg]

%% Gauss-Laguerre Quadrature (accurate up to 2n-1)
% Represents:
% \int_0^{\inf} f(x) exp(-x) dx = \sum_i wi f(xi)

Nq = 5;
[xi, wi] = LAGWT(Nq);

[pl, intg] = PLAGU(0, xi);

[wi'*(pl.^2) intg]

%% Gauss-Jacobi Quadrature
% Represents (defaults to a=0, b=1
% \int_{a}^{b} f(x) (b-x)^(bt)*(x-a)^(alf) dx = \sum_i wi f(xi)
% \int_{0}^{1} f(x) x^(alf)*(1-x)^(bt) dx = \sum_i wi f(xi)
alf = 3.5;
bt = 2.5;

Nq = 5;
[xi, wi] = GJWT(Nq, alf, bt, -1,1);

[pl, intg] = PJACO(0, alf,bt, xi, -1,1);

[wi'*(pl.^2) intg]
