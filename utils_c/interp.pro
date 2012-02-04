
N = 1.e6

xfull = dindgen(N)
yfull = sin(xfull/3.1415926535897931d0/2/20)+0.2
xint = dindgen(11)/10.*N

t0 = systime(/seconds)

yint_0 = interpol(yfull, xfull, xint)
t1 = systime(/seconds)

print, t1-t0

print, sin(xint/3.1415926535897931d0/2/20)+0.2-yint_0