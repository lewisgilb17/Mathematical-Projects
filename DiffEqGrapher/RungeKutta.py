import math

def RK4_2(x, y, z, dx, dydx, dzdx):
    
    k1 = dx*dydx(x, y, z)
    h1 = dx*dzdx(x, y, z)
    k2 = dx*dydx(x+dx/2., y+k1/2., z+h1/2.)
    h2 = dx*dzdx(x+dx/2., y+k1/2., z+h1/2.)
    k3 = dx*dydx(x+dx/2., y+k2/2., z+h2/2.)
    h3 = dx*dzdx(x+dx/2., y+k2/2., z+h2/2.)
    k4 = dx*dydx(x+dx, y+k3, z+h3)
    h4 = dx*dzdx(x+dx, y+k3, z+h3)

    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    z = z + 1./6.*(h1+2*h2+2*h3+h4)    
    return y, z

def RK4_1(x, y, dx, dydx):
    
    # Calculate slopes
    k1 = dx*dydx(x, y)
    k2 = dx*dydx(x+dx/2., y+k1/2.)
    k3 = dx*dydx(x+dx/2., y+k2/2.)
    k4 = dx*dydx(x+dx, y+k3)
    
    # Calculate new x and y
    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    
    return y

def RK4_3(x, y, z, w, dx, dydx, dzdx, dwdx):
    
    k1 = dx*dydx(x, y, z, w)
    h1 = dx*dzdx(x, y, z, w)
    n1 = dx*dwdx(x, y, z, w)

    k2 = dx*dydx(x+dx/2., y+k1/2., z+h1/2., w+n1/2.)
    h2 = dx*dzdx(x+dx/2., y+k1/2., z+h1/2., w+n1/2.)
    n2 = dx*dwdx(x+dx/2., y+k1/2., z+h1/2., w+n1/2.)

    k3 = dx*dydx(x+dx/2., y+k2/2., z+h2/2., w+n2/2.)
    h3 = dx*dzdx(x+dx/2., y+k2/2., z+h2/2., w+n2/2.)
    n3 = dx*dwdx(x+dx/2., y+k2/2., z+h2/2., w+n2/2.)

    k4 = dx*dydx(x+dx, y+k3, z+h3, w+n3)
    h4 = dx*dzdx(x+dx, y+k3, z+h3, w+n3)
    n4 = dx*dwdx(x+dx, y+k3, z+h3, w+n3)

    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    z = z + 1./6.*(h1+2*h2+2*h3+h4)   
    w = w + 1./6.*(n1+2*n2+2*n3+n4) 
    return y, z, w

def RK4_4(x, y, z, w, v, dx, dydx, dzdx, dwdx, dvdx):
    
    k1 = dx*dydx(x, y, z, w, v)
    h1 = dx*dzdx(x, y, z, w, v)
    n1 = dx*dwdx(x, y, z, w, v)
    l1 = dx*dvdx(x, y, z, w, v)

    k2 = dx*dydx(x+dx/2., y+k1/2., z+h1/2., w+n1/2., v+l1/2.)
    h2 = dx*dzdx(x+dx/2., y+k1/2., z+h1/2., w+n1/2., v+l1/2.)
    n2 = dx*dwdx(x+dx/2., y+k1/2., z+h1/2., w+n1/2., v+l1/2.)
    l2 = dx*dvdx(x+dx/2., y+k1/2., z+h1/2., w+n1/2., v+l1/2.)

    k3 = dx*dydx(x+dx/2., y+k2/2., z+h2/2., w+n2/2., v+l2/2.)
    h3 = dx*dzdx(x+dx/2., y+k2/2., z+h2/2., w+n2/2., v+l2/2.)
    n3 = dx*dwdx(x+dx/2., y+k2/2., z+h2/2., w+n2/2., v+l2/2.)
    l3 = dx*dvdx(x+dx/2., y+k2/2., z+h2/2., w+n2/2., v+l2/2.)

    k4 = dx*dydx(x+dx, y+k3, z+h3, w+n3, v+l3)
    h4 = dx*dzdx(x+dx, y+k3, z+h3, w+n3, v+l3)
    n4 = dx*dwdx(x+dx, y+k3, z+h3, w+n3, v+l3)
    l4 = dx*dvdx(x+dx, y+k3, z+h3, w+n3, v+l3)

    y = y + 1./6.*(k1+2*k2+2*k3+k4)
    z = z + 1./6.*(h1+2*h2+2*h3+h4)   
    w = w + 1./6.*(n1+2*n2+2*n3+n4) 
    v = v + 1./6.*(l1+2*l2+2*l3+l4) 

    #return graph of derivative last
    derivative1 = k1 / dx
    derivative2 = h1 / dx
    return y, z, w, v, derivative1, derivative2