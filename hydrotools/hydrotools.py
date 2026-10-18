"""Main module."""

from math import atan2,cos, sin, degrees

# Determine if a point is inside a given polygon or not
# Polygon is a list of (x,y) pairs. This fuction
# returns True or False.  The algorithm is called
# "Ray Casting Method".
def point_in_poly(x,y,poly):
    if (x, y) in poly: return True
    for I in range(len(poly)):
        p1 = None
        p2 = None
        if I==0:
            p1 = poly[0] 
            p2 = poly[1]
        else:
            p1 = poly[I-1]
            p2 = poly[I]
        if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
            return True
        
        n = len(poly)
        inside = False

        p1x,p1y = poly[0]
        for i in range(n+1):
            p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y

    return inside

# Calculate line direction
def bearing(x1, y1, x2, y2):
    angle = atan2(cos(y1)*sin(y2)-sin(y1)*cos(y2)*cos(x2-x1),sin(x2-x1)*cos(y2))
    return (degrees(angle)+360)%360