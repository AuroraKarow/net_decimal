# note

## reciprocal

$$ x_{n+1}=x_n\frac{f(x_n)}{f'(x_n)} $$
$$ a=g(x)=\frac{1}{x}(x\in\R) $$
$$ f(x)=g(x)-a=\frac{1}{x}-a=0 $$
$$ x_{n+1} = x_n(2-x_na) $$

## exponent

$$ e^x=1+x+\frac{x^2}{2!}+\frac{x^3}{3!}+\dots+\frac{x^n}{n!}(x\in\R) $$
$$ c_0=1,y_0=c,s_0=c,x_0=x $$
$$ \begin{align*}
    y_{n+1}&=y_n+\frac{x_n}{s_n} \\
    x_{n+1}&=x_nx \\
    c_{n+1}&=c_n+1 \\
    s_{n+1}&=s_nc_{n+1}
\end{align*} $$

## logarithm

$$ \ln x = 2\left[\left(\frac{x-1}{x+1}\right)+\frac{1}{3}\left(\frac{x-1}{x+1}\right)^3+\frac{1}{5}\left(\frac{x-1}{x+1}\right)^5+\dots+\frac{1}{2n+1}\left(\frac{x-1}{x+1}\right)^{2n+1}\right](x>0) $$
$$ b_0=1,s_0=\frac{x-1}{x+1},y_0=0,m=s_0^2 $$
$$ \begin{align*}
    y_{n+1}&=y_n+\frac{s_n}{b_n} \\
    s_{n+1}&=s_nm \\
    b_{n+1}&=b_n+2 \\
    \ln x&=2y_n
\end{align*} $$
$$ \ln(x+1)=x-\frac{x^2}{2}+\frac{x^3}{3}-\frac{x^4}{4}+\dots+(-1)^n\frac{x^{n+1}}{n+1},(-1<x\le1) $$
$$ b_0=1,y_0=0,c_0=b,x_0=x $$
$$ \begin{align*}
    y_{n+1}&=y_n+c_0\frac{x_n}{b_n} \\
    c_{n+1}&=-c_n \\
    b_{n+1}&=b_n+1 \\
    x_{n+1}&=x_{n}x
\end{align*} $$

## sine

$$ \sin x=x-\frac{x^3}{3!}+\frac{x^5}{5!}-\frac{x^3}{3!}+\dots+(-1)^n\frac{x^{2n+1}}{(2n+1)!}(x\in\R) $$
$$ y_0=0,x_0=x,b_0=1,s_0=b_0,c_0=b_0 $$
$$ \begin{align*}
    y_{n+1}&=y_n+c_n\frac{x_n}{s_n} \\
    s_{n+1}&=s_n(b_0+1)(b_0+2) \\
    b_{n+1}&=b_n+2 \\
    x_{n+1}&=x_nx^2 \\
    c_{n+1}&=-c_n
\end{align*} $$
$$ \begin{align*}
    \sin(x+2k\pi)&=\sin x(k\in\Z) \\
    \sin(\pi+x)&=-\sin x \\
    \sin(\pi-x)&=\sin x \\
    \sin(\frac{\pi}{2}+x)&=\cos x \\
    \sin(\frac{\pi}{2}-x)&=\cos x
\end{align*} $$

## cosine

$$ \cos x=1-\frac{x^2}{2!}+\frac{x^4}{4!}-\frac{x^6}{6!}+\dots+(-1)^n\frac{x^{2n}}{(2n)!}(x\in\R) $$
$$ y_0=0,x_0=1,b_0=y_0,s_0=x_0,c_0=x_0 $$
$$ \begin{align*}
    y_{n+1}&=y_n+c_n\frac{x_n}{s_n} \\
    s_{n+1}&=s_n(b_0+1)(b_0+2) \\
    b_{n+1}&=b_n+2 \\
    x_{n+1}&=x_nx^2 \\
    c_{n+1}&=-c_n \\
\end{align*} $$
$$ \begin{align*}
    \cos(x+2k\pi)&=\cos x(k\in\Z) \\
    \cos(\pi+x)&=\cos x \\
    \cos(\pi-x)&=-\cos x \\
    \cos(\frac{\pi}{2}+x)&=-\sin x \\
    \cos(\frac{\pi}{2}-x)&=\sin x
\end{align*} $$

## $\pi$

$$ \pi=\sum_{i=0}^n\frac{1}{16^i}\left(\frac{4}{8i+1}-\frac{2}{8i+4}-\frac{1}{8i+5}-\frac{1}{8i+6}\right) $$
$$ p=8,q=\frac{1}{16},c_0=1,d_0=0,a_0=d_0 $$
$$ \begin{align*}
    a_{n+1}&=a_n+c_n\left(\frac{4}{d_n+1}-\frac{2}{d_n+4}-\frac{1}{d_n+5}-\frac{1}{d_n+6}\right) \\
    c_{n+1}&=c_nq \\
    d_{n+1}&=d_n+p \\
\end{align*} $$

## power

$$ f(x)=x^a(x,a\in\R) $$
$$ f(x)=\begin{cases}
    0&x=0\\
    e^{a\ln x}&x>0\\
    |x|^ae^{i(2k+1)\pi a}&x<0,k\in\Z
\end{cases} $$
$$ e^{ix}=\cos x+i\sin x $$
$$ |x|^ae^{i(2k+1)\pi a}=|x|^a(\cos T+i\sin T) $$
$$ T=(2k+1)\pi a,(k\in\Z) $$
