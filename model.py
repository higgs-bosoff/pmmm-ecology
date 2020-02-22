import numpy as np

class MWModel(object):

    def __init__(self, r=1, t0=1, b=0.2, w=0.5, eta=0.25):
        
        self.r = r 
        self.t0 = t0
        self.b = b
        self.w = w 
        self.eta = eta

    def set_population(self, M, W):

        self.M = M 
        self.W = W

    def decayRate(self):
        return 1.0/(self.t0+self.eta*self.b*self.w*self.W)

    def eqWpop(self):
        return ((self.t0**2+4*self.eta)-self.t0)/(2*self.eta*self.b*self.w)

    def eqTotDecayRate(self):
        W = self.eqWpop()

        return 1/(self.t0+self.eta*self.b*self.w*W)+self.b*(1-self.w)

    def step(self, dt):

        dM = self.r-self.decayRate()*self.M-self.b*(1-self.w)*self.M*self.W
        dW = self.decayRate()*self.M-self.b*self.w*self.M*self.W

        self.M += dM*dt
        self.W += dW*dt