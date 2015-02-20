library(manipulate)
library(flowCore)

fcs.data <- read.FCS('~/thor/FCS/exvivo-IL2/CB00366X_01U_2012-11-07.fcs')

load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00366X_2012-11-07.RData' )
load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB01510Q_2012-11-21.RData')

#APC-CD25
x <- fcs.data@exprs[,9]
plot(density(x))

#w <- (m-log10(t/abs(r)))/2
# r is the most negative value to be included in the display

# @param w is the linearization width in asymptotic decades. w should be > 0 and determines the slope of transformation at zero.
# w can be estimated using the equation w=(m-log10(t/abs(r)))/2, where r is the most negative value to be included in the display
# @param t	Top of the scale data value, e.g, 10000 for common 4 decade data or 262144 for a 18 bit data range. t should be greater than zero
# @param m is the full width of the transformed display in asymptotic decades. m should be greater than zero
# @param a	Additional negative range to be included in the display in asymptotic decades. Positive values of the argument brings additional negative input values into the transformed display viewing area. Default value is zero corresponding to a Standard logicle function.
trans <- function(w=w, t=t, m = m, a = a) logicleTransform(w=w, t = t, m = m, a = a)

T*10**-(m-w-a) * (10**(x-w-a) - r**2 * 10**-((x-w-a)/r) + r**2 - 1)
-T*10**-(m-w-a) * (10**(w+a-x) - r**2 * 10**-((w+a-x)/r) + r**2 - 1)


load('~/thor/CB01510Q_1000U_2012-11-21.RData')
x2 <- fcs.data[,'PSTAT5']

x <- fcs.data[,'PSTAT5.4']
manipulate( plot(density(trans(w,t=262144,m,a)(x))), w=slider(0,10,initial=.5, step=.1), m=slider(0,10,initial=4.5, step=.5), a=slider(-5,5,initial=0) )

manipulate( plot(density(trans(w,t=262144,m,a)(x2))), w=slider(0,10,initial=.5, step=.1), m=slider(0,10,initial=4.5, step=.5), a=slider(-5,5,initial=0) )


f <- function(w,t,m,a) {
  curve(trans(w,t=262144,m,a)(x),from=-10**4,to=10**4)
  curve(trans(w=0,t=262144,m=4.5,a=0)(x), from=-10**4, to=10**4,lwd=2,add=TRUE)
  abline(v=0)
  abline(h=0)
}
manipulate(f(w,t,m,a), w=slider(0,10,initial=.5, step=.1), m=slider(0,10,initial=4.5, step=.5), a=slider(-5,5,initial=0) )
#like asinh but with power of 10 instead of e
#curve(asinh(x),from=-10**4,to=10**4)
