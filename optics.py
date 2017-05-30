fc = 2
wd = 1080/fc
ht = 1080/fc
ctrx = wd/2.   # center
ctry = ht/2.
### circle medium
mediumc = PVector(ctrx,ctry)  # circular medium center
mediumr = ctrx/2       # circular medium radius
### plane medium
mediump = PVector(0,ht) #direction
mediumd = PVector(ctrx,0) #displacement
mediuml = wd/5. # width
###
nwl = 3  # n of wavelengths
npts = 128  # n of particles per wavelength
concentric = False   # concentric / planar photons
circularMedium = True   # circular / planar medium
edges = False
### CONCENTRIC PHOTONS
px = ctrx*.5           # coords for concentric photons
py = ctry
anglerange = (PI/4,PI)  # angle range for concentric photons
### LINEARLY PHOTONS
# th1t = 41/2.
# b = mediumr*1.33*sin(th1t)
xif = (ctrx*.51,ctrx*.99)  # initial/final values for x por linear position
def lineeq(x):   # line eq that returns y given x depending on endpoints
    p1 = PVector(0,0)
    p2 = PVector(ctrx*2,0)
    y = (p2.y-p1.y)/(p2.x-p1.x)*(x-p1.x) + p1.y
    return y
dir = PI/2   # velocity direction of linearly arranged photons
### etc
strokerange = (2./fc,40./fc)  # range of stroke width
swidth = 3  # stroke width of photons
spd = 1     # initial speed of photons
Ilimit = 0.0   # minimal photon intensity before it is destroyed
fadeamt = 0.1

def settings():
    size(wd,ht)
    
def setup():
    global spectrum
    colorMode(HSB,1.0)
    ellipseMode(RADIUS)
    background(0)
    noFill()
    spectrum = []
    for i in xrange(nwl):
        off = 1
        if nwl == 1: off = 0
        wl = map(i,0,nwl-off,.4,.7)
        s = rays(npts,wl,concentric,circularMedium)
        spectrum.append(s)

def draw():
    # background(0)
    blendMode(SUBTRACT)
    noStroke()
    fill(fadeamt)
    rect(0,0,width,height)
    blendMode(ADD)
    noFill()
    stroke(0,fadeamt)
    strokeWeight(2)
    ellipse(mediumc.x,mediumc.y,mediumr,mediumr)
    # stroke(1)
    # strokeWeight(1)
    # line(mediumd.x-mediuml,0,mediumd.x-mediuml,ht)
    # line(mediumd.x+mediuml,0,mediumd.x+mediuml,ht)
    for i in xrange(len(spectrum)):
        spectrum[i].step()

###
class rays:
    def __init__(self,n,wl,concentric,circle):
        self.n = n   # number of ray groups
        self.wavelength = wl  # wavelength of group
        self.nref = nlambda(wl)  # refraction index depending on wl
        h = map(wl,.4,.7,.75,0)  # wavlenght mapped to hue
        self.c = color(h,1,1)
        self.p = []
        for i in xrange(n):
            # concentric photons
            if concentric:
                theta = map(i,0,n,anglerange[0],anglerange[1])
                vel = PVector.fromAngle(theta)
                vel.setMag(spd)
                P = particle(PVector(px,py),vel,self.nref,circle)
            # linearly arranged photons
            else:
                off = 1
                if n == 1: off = 0
                x = map(i,0,n-off,xif[0],xif[1])
                y = lineeq(x)
                pos = PVector(x,y)
                vel = PVector.fromAngle(dir)
                vel.setMag(spd)
                P = particle(pos,vel,self.nref,circle)
            self.p.append(P)
    
    def update(self):
        newparticles = []  #photons created at boundary
        deadparticles = [] #photons that need to be delet
        for i in xrange(len(self.p)):
            # updates the current photon,
            # determines if a new one needs to be created
            # of if an existing one needs to be destroyed
            new,kill = self.p[i].update()
            if kill:
                deadparticles.append(i)
            if new != None:
                newparticles.append(new)
        # remove photons
        if len(deadparticles) != 0:
            for i in xrange(len(deadparticles)-1,-1,-1):
                self.p.pop(deadparticles[i])
        # add photons
        if len(newparticles) != 0:
            self.p.extend(newparticles)

    def show(self):
        strokeWeight(swidth)
        for i in xrange(len(self.p)):
            # value depending on photon intensity
            sw = map(self.p[i].I,Ilimit,1,strokerange[0],strokerange[1])
            alp = map(self.p[i].I,Ilimit,1,.1,1)
            strokeWeight(sw)
            stroke(hue(self.c),1,1,alp)
            line(self.p[i].prevpos.x,self.p[i].prevpos.y,
                 self.p[i].pos.x,self.p[i].pos.y)
    
    def step(self):
        self.update()
        self.show()
###
###

class particle:
    def __init__(self,pos,vel,nref,circle):
        self.pos = pos  # current position
        self.prevpos = pos # previous position to draw line
        self.vel = vel  # velocity vector
        self.nref = nref # refraction index
        self.I = 1  # photon "intensity"
        self.circle = circle
        # calculate refraction index based on position
        ### circle
        if self.circle:
            dif = PVector.sub(self.pos,mediumc)
            if dif.mag() < mediumr: self.n1 = self.nref
            else: self.n1 = 1
        ### plane
        else:
            if (mediumd.x-mediuml) < self.pos.x < (mediumd.x+mediuml):
                self.n1 = self.nref
            else: self.n1 = 1    
    
    def update(self):
        new = None
        kill = False
        self.prevpos = self.pos.copy()
        # calculate refraction index based on position
        ### circle
        if self.circle:
            dif = PVector.sub(self.pos,mediumc)
            if dif.mag() < mediumr: n2 = self.nref
            else: n2 = 1
        ### plane
        else:
            if (mediumd.x-mediuml) < self.pos.x < (mediumd.x+mediuml):
                n2 = self.nref
            else: n2 = 1
        # refraction index ratio
        nratio = self.n1/n2
        if nratio != 1:
            ### circle boundary
            if self.circle:
                nvec = PVector.sub(self.pos,mediumc)
            ### plane boundary
            else:
                pv = PVector.sub(self.pos,mediumd)
                b = mediump.copy()
                b.normalize()
                aa1 = PVector.dot(pv,b)
                a1 = PVector.mult(b,aa1)
                nvec = PVector.sub(pv,a1)
            nvec.normalize()
            if nratio > 1:
                nvec.mult(-1)
            #
            ivec = self.vel.copy()
            ivec.normalize() # normalized incident vector
            cosi = -PVector.dot(nvec,ivec)
            sin2t = nratio**2 * (1 - cosi**2)
            ### REFLECTANCE
            rft = 1   # reflectance
            if sin2t <= 1:
                r0 = ((self.n1-n2)/(self.n1+n2))**2
                cosx = cosi
                if nratio > 1:
                    cosx = (1 - sin2t)**.5
                rft = r0 + (1-r0)*(1-cosx)**5
            ### REFLECT
            if rft >= Ilimit:
                # calculate/create reflected vector
                cosn2 = PVector.mult(nvec,2*cosi)
                rflvec = PVector.add(ivec,cosn2)
                newvec = rflvec.copy()
                newvec.setMag(spd)
                newpos = PVector(self.pos.x+rflvec.x,self.pos.y+rflvec.y)
                new = particle(newpos,newvec,self.nref,self.circle)
                new.I *= rft
                # new.pos.add(new.vel)
            ### REFRACT
            if sin2t <= 1:
                # calculate/delete refracted vector
                in1n2 = PVector.mult(ivec,nratio)
                nn1n2 = PVector.mult(nvec,nratio*cosi - (1 - sin2t)**.5)
                rfrvec = PVector.add(in1n2,nn1n2)
                ### apply
                vmag = self.vel.mag()
                self.vel = rfrvec.copy()
                self.vel.setMag(myround(vmag*nratio,.0001))
                self.I *= 1 - rft
                if self.I < Ilimit: kill = True
            else:
                kill = True
        self.n1 = n2 # update new refraction index
        self.pos.add(self.vel) # update position
        out = self.checkBoundary(edges) # check for oob conditions
        if out: kill = True
        return new,kill

    def checkBoundary(self,boundary=True):
        out = False
        if not(0 <= self.pos.x < width):
            if boundary:
                self.vel.x *= -1
            else: out = True
        if not(0 <= self.pos.y < height):
            if boundary:
                self.vel.y *= -1
            else: out = True
        return out

###

def nlambda(wl):
    ### BK7 GLASS
    # B1 = 1.0396
    # B2 = 0.23179
    # B3 = 1.0105
    # C1 = .0060007
    # C2 = .020018
    # C3 = 103.56
    ### WATER 19C
    B1 = .56725
    B2 = .17366
    B3 = .021215
    C1 = .0050856
    C2 = .018149
    C3 = .026172
    n2 = 1 + B1*wl**2/(wl**2-C1) + B2*wl**2/(wl**2-C2) + B3*wl**2/(wl**2-C3)
    return n2**.5

def myround(x, base):
    return base * round(float(x)/base)

def keyPressed():
    try:
        if key.lower() == 'f': print frameCount
    except:
        pass
