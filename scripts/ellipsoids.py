"""
Defines classes and functions to study couples of ellipsoids.
"""

from sage.matrix.constructor import matrix as Matrix
from sage.calculus.var import var
from sage.calculus.functional import expand
from sage.functions.other import sqrt

import numpy as np

from numbers import Number

##################
### ELLIPSOIDS ###
##################

class Ellipsoid:
    """
    Description of a single ellipsoid.
    """

    def __init__(self, x, y, z, q0, q1, q2, q3, R1, R2, R3):
        """
        Initiate position, quaternion, and semi-axes of ellipsoid.
        """

        self.pos = Matrix([[x], [y], [z]])
        self.quat = Quaternion(q0, q1, q2, q3)
        self.R = np.array([R1, R2, R3])

    def rotation_matrix(self):
        """
        Returns rotational matrix associated to the quaternion of the ellipsoid.
        """

        q0, q1, q2, q3 = self.quat.coords
        Q = Matrix(
            [[1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
            [2*(q1*q2 + q0*q3), 1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
            [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), 1 - 2*(q1**2 + q2**2)]])
        return Q

    def red_belonging_matrix(self):
        """
        Returns reduced belonging matrix associated to ellipsoid.
        """

        Q = self.rotation_matrix()
        R1, R2, R3 = self.R
        diag = Matrix([[1/(R1**2), 0, 0], [0, 1/(R2**2), 0], [0, 0, 1/(R3**2)]])
        return Q*diag*Q.transpose()

    def belfunc(self):
        """
        Returns belonging function of ellipsoid.
        """

        X = Matrix([[var('x')], [var('y')], [var('z')]])
        v = self.pos
        B = self.red_belonging_matrix()
        return (expand((X - v).transpose()*B*(X - v)) - 1)[0, 0]

    def rescale_axes(self, factor):
        """
        Rescale axes by factor factor.
        """

        self.R *= factor
        return self

    def rotate(self, axis, angle):
        """
        Add rotation of angle angle along axis described by array-like axis of
        length 3 to the quaternion.
        """

        self.quat.rotate(axis, angle)
        return self

class EllipsoidCouple:
    """
    Description of a couple of ellipsoid.
    """

    def __init__(self, ellipsoid1, ellipsoid2):
        """
        Initiate couple of ellipsods.
        """

        self.ellipsoid1 = ellipsoid1
        self.ellipsoid2 = ellipsoid2

    def mu(self, epsilon=1e-12):
        """
        Returns rescaling factor associated to the couple with precision epsilon
        in the search of the root of the contact function.
        """

        eta = var('eta')

        Y = (eta*self.ellipsoid2.red_belonging_matrix().inverse() +
            (1 - eta)*self.ellipsoid1.red_belonging_matrix().inverse())

        p = (eta*(1 - eta)*
            ((self.ellipsoid2.pos - self.ellipsoid1.pos).transpose())*
            (Y.adjugate())*(self.ellipsoid2.pos - self.ellipsoid1.pos))
        q = Y.determinant()

        h = expand(p.derivative(eta)*q - p*q.derivative(eta))[0, 0]
        F = expand(p/q)[0, 0]

        root, iter = halleyMethod(
            h, eta,
            self.ellipsoid1.R[0]/(self.ellipsoid1.R[0] + self.ellipsoid2.R[0]),
            epsilon)
        return sqrt(F(eta=root))

#############
### MATHS ###
#############

def halleyMethod(f, x, x0, epsilon=1e-12):
    """
    Returns estimate of a root of the function f of x estimated with Halley's
    method with initial guess x0 and precision epsilon, as well as the number
    of iterations of the method.
    """

    fp = f.derivative(x)    # first derivative of f with respect to x
    fpp = fp.derivative(x)  # second derivativate of f with respect to x

    iterations = 0              # number of iterations
    while abs(f(x=x0)) > epsilon: # Hayley's method iteration
        iterations += 1
        x0 -= (f(x=x0)/fp(x=x0))/(
            1 - (f(x=x0)/fp(x=x0))*(fpp(x=x0)/(2*fp(x=x0))))

    return x0, iterations

class Quaternion:
    """
    Use quaternions.
    """

    def __init__(self, q0=1, q1=0, q2=0, q3=0):
        """
        Initiates coordinates of quaternion.
        """

        self.coords = np.array([q0, q1, q2, q3], dtype=float)

    def __mul__(self, other):
        """
        Defines quaternion multiplication.
        """

        coords = np.empty(4, dtype=float)   # new coordinates of the quaternion

        if isinstance(other, Quaternion):   # multiplication with an other quaternion
            coords[1:], coords[0] = (
                (np.cross(self.coords[1:], other.coords[1:])
                    + self.coords[0]*other.coords[1:]
                    + other.coords[0]*self.coords[1:]),
                (self.coords[0]*other.coords[0]
                    - np.dot(self.coords[1:], other.coords[1:])))
        elif isinstance(other, Number):     # multiplication with a scalar
            coords = other*self.coords
        else:
            raise TypeError(f"Cannot multiply by {type(other).__name__}")

        return Quaternion(*coords)

    def rotate(self, axis, angle):
        """
        Add rotation of angle angle along axis described by array-like axis of
        length 3 to the quaternion.
        """

        self.coords = (
            Quaternion(np.cos(angle/2), *(np.array(axis)*np.sin(angle/2)))
            *self).coords
        return self

    def normalise(self):
        """
        Normalise coordinates.
        """

        self.coords = self.coords/np.sqrt(self.coords**2)

    def __repr__(self):
        """
        Prints quaternion coordinates.
        """

        return str(self.coords.tolist())
