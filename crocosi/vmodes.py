# -*- coding: utf-8 -*-
"""
Vertical mode projection class

"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
from IPython import embed
from utils import tvs_to_s_fast, rhoz, volint
import os


class VModes(object):
    """
    Vertical modes class
    """
    def __init__(self, zc, zf, N2f, nmodes, g=9.81, free_surf=True, sigma=0.1):
        """
        Constructor, computes the vertical modes by calling _computemodes()
        zc are the cell centres, zf are the cell faces, both in ascending order
        N2f is N^2(zf), nmodes is number of modes
        """
        self.zc, self.zf, self.N2f, self.nmodes, self.g = zc, zf, N2f, nmodes, g
        self.N2c = np.interp(self.zc,self.zf,self.N2f) # N2 on c points
        self.H   = np.abs(self.zf[0])
        self.free_surf = free_surf
        self.sigma=sigma
        self._computemodes()
        
    def _computemodes(self):
        # Precompute a few quantities
        dzc=np.diff(self.zc)
        dzf=np.diff(self.zf)
        Np=len(self.zf)
        N20=self.N2f[-1]  # N^2(z=0)
        
        # Build Dz, C-to-F grid        
        v1=-1.0/np.append(dzc, 1)
        v2= 1.0/np.append(1, dzc)
        v12=np.append([v1],[v2],axis=0)
        Dz=sp.spdiags(v12,[-1, 0],Np,Np-1,format="lil")
        # Adjust matrix for BCs
        Dz[0,:]=0
        Dz[-1,:]=0
        if self.free_surf:
            Dz[-1,-1]=np.divide(-N20, self.g + N20*(self.zf[-1] - self.zc[-1]))
        
        # Build Dz2, F-to-C grid
        v1=-1.0/np.append(dzf,1)
        v2= 1.0/np.append(1,dzf)
        v12=np.append([v1],[v2],axis=0)
        Dz2=sp.spdiags(v12,[0, 1],Np-1,Np,format="lil")
        
        # Construct A, solve eigenvalue problem
        iN2=sp.spdiags(1.0/self.N2f,0,Np,Np)
        A=-Dz2*iN2*Dz
        #ev,ef = la.eigs(A.tocsc(),self.nmodes+1,sigma=0.1)
        ev,ef = la.eigs(A.tocsc(),self.nmodes+1,sigma=self.sigma)
        ev,ef = np.real(ev), np.real(ef)
        
        # Convert eigvenvalues to c_e, sort appropriately
        c=1.0/np.sqrt(np.real(ev))
        ii=(-c).argsort()
        self.c=c[ii]               # c_e value
        self.phic=ef[:,ii]         # phi at cell centres

        # Normalize and set surface value positive
        for mi in range(self.nmodes+1):
            fn=self.phic[:,mi]        # current phi
            s=np.sign(fn[-1])         # sign at surface
            if s==0:
                s=1;
            tmp = np.sum((fn**2)*dzf)/self.H
            self.phic[:,mi] = s*fn/np.sqrt(tmp) # (1/H)*\int_{-H}^{0} \phi_m(z)^2 dz = 1

        # Verify number of zero crossings matches mode number
        for mi in range(self.nmodes+1):
            nzc=int(np.sum(np.abs(np.diff(np.sign(self.phic[:,mi]))))/2.0)
            if mi != nzc:
                print("problem! starting embed() in vmodes._computemodes()")
                embed()
                exit()

        # dphi at cell faces
        self.dphif=Dz*self.phic

        # Linear interpolate dphi at cell faces to find dphi at cell centres
        self.dphic = np.zeros(self.phic.shape)
        for m in range(self.nmodes+1):
            self.dphic[:,m] = np.interp(self.zc, self.zf, self.dphif[:,m])

        self.orthogonality()

    def orthogonality(self):
        """
        Check orthogonality
        """
        dzf    = np.diff(self.zf)
        N40    = self.N2f[-1]**2  # = N^4(z=0)
        phio_e = np.eye(self.nmodes+1)                    # analytic expectation
        phio   = np.zeros([self.nmodes+1,self.nmodes+1])  # computed orthogonality
        dphio_e= np.zeros([self.nmodes+1,self.nmodes+1])  # analytic expectation
        dphio  = np.zeros([self.nmodes+1,self.nmodes+1])  # computed orthogonality
        dzfoN2 = dzf/self.N2c
        for n in range(self.nmodes+1):
            cn, phicn, dphicn, dphifn = self.getmode(n)
            for m in range(self.nmodes+1):
                cm, phicm, dphicm, dphifm = self.getmode(m)
                # phi orthogonality
                phio[n,m]    = np.sum(phicn*phicm*dzf)/self.H
                # dphi orthogonality
                dphio[n,m]   = np.sum(dphicn*dphicm*dzfoN2)/self.H
                # dphi expectation
                dphio_e[n,m] = phio_e[n,m]/(cn*cn) -self.g*dphifn[-1]*dphifm[-1]/(N40*self.H)

        if False:
            # Also evaluate alternate approach to find dphic
            dphic2 = np.zeros(self.phic.shape)
            for m in range(self.nmodes+1):  # dphic2: use utils.rhoz
                dphic2[:,m] = rhoz(self.phic[:,m][:,np.newaxis,np.newaxis],self.zc)[:,0,0]
            dphic3       = np.copy(dphic2)   # dphic3: same as dphic2, but
            dphic3[0 ,:] = self.dphic[0 ,:]  #         use linear interpolated
            dphic3[-1,:] = self.dphic[-1,:]  #         endpoints from dphic

            dphio2 = np.zeros([self.nmodes+1,self.nmodes+1])  # computed orthogonality (dphic2)
            dphio3 = np.zeros([self.nmodes+1,self.nmodes+1])  # computed orthogonality (dphic3)
            for n in range(self.nmodes+1):
                for m in range(self.nmodes+1):
                    # alternate dphi cases
                    dphio2[n,m]=np.sum(dphic2[:,n]*dphic2[:,m]*dzfoN2)/self.H
                    dphio3[n,m]=np.sum(dphic3[:,n]*dphic3[:,m]*dzfoN2)/self.H
            e1 = np.abs((dphio_e-dphio ))
            e2 = np.abs((dphio_e-dphio2))
            e3 = np.abs((dphio_e-dphio3))
            embed()

    def getmode(self,m):
        """
        Extract a mode, return ce,phic,dphic,dphif tuple
        m=0 returns barotropic mode
        m=1,2,3,... returns baroclinic modes
        """
        return self.c[m], self.phic[:,m], self.dphic[:,m], self.dphif[:,m]

    def tripleproducts(self):
        """
        Computes the triple products
        """
        siz=[self.nmodes+1, self.nmodes+1, self.nmodes+1] # output array size
        self.alpha, self.beta  = np.zeros(siz), np.zeros(siz)
        self.gamma, self.delta = np.zeros(siz), np.zeros(siz)
        self.lambd, self.pi    = np.zeros(siz), np.zeros(siz) # can't use lambda
        dzf = np.diff(self.zf)
        for n in range(self.nmodes+1):
            cn, phicn, dphicn, dphifn = self.getmode(n)
            for m in range(self.nmodes+1):
                cm, phicm, dphicm, dphifm = self.getmode(m)
                ddphicm=np.diff(dphifm)/dzf
                for l in range(self.nmodes+1):
                    cl, phicl, dphicl, dphifl = self.getmode(l)
                    self.alpha[n,m,l] = np.sum(phicn*phicm*phicl*dzf)/self.H
                    self.beta [n,m,l] = np.sum(phicn*dphicm*dphicl*dzf/self.N2c)/self.H
                    self.gamma[n,m,l] = np.sum(dphicn*dphicm*phicl*dzf*cn*cn/self.N2c)/self.H
                    self.delta[n,m,l] = np.sum(dphicn*ddphicm*dphicl*dzf*cn*cn/(self.N2c**2))/self.H
                    self.lambd[n,l,m] = self.gamma[n,m,l]
                    self.pi   [n,l,m] = self.delta[n,m,l]
        return self.alpha, self.beta, self.gamma, self.delta, self.lambd, self.pi


class N2Profile(object):
    """
    Class to compute and hold the N2 profile
    """
    def __init__(self, case, ti, regs, g=9.81, denfile=None):
        """
        Handles the N2 profiles
          1) Loads the right average density file
          2) Computes spatial average density profiles
          3) Computes N2 profiles
        """
        from utils import horizmean
        self.zf = case.getZ('w')        # Nominal z grid at w points, 1D array
        self.zc = case.getZ('t')        # Nominal z grid at t points, 1D array
        yT = case.getXY('T')[1]; yT=yT[:,0]   # y values of the V points

        # Load density for constructing a stratification profile
        if denfile is None:
            # Use the 1-day averaged density rho_a from model outputs
            rin = case.his.variables['T_a'][ti,:,:,:]
            ssh = case.his.variables['ssh_a'][ti,:,:]
            self.r = tvs_to_s_fast(rin,ssh,case)
        else:
            # Load the time-averaged and interpolated density from denfile
            npzfile = np.load(denfile)
            self.r = npzfile['r']

        # Create empty dict()s to store N2 profiles
        self.rr, self.N2 = dict(), dict()
        for ri,val in enumerate(regs):
            if type(val) is tuple:       # Compute the bounded domains
                ymin,ymax = val
                ir = np.where( (yT >  ymin) & (yT < ymax))[0]  # indices for rho
                self.rr[val] = np.mean(np.mean(self.r[:,ir,:], axis=1),axis=1)
            # Compute the presets
            elif val == "south":
                self.rr[val] = self.r[:,1,1]
            elif val == "north":
                self.rr[val] = self.r[:,-2,1]
            elif val == "mean":
                self.rr[val] = horizmean(self.r,case.hgrid)
            elif val == "middle-third":
                self.rr[val] = np.mean(np.mean(self.r[:,(case.Mm/2-case.Mm/6):(case.Mm/2+case.Mm/6),1:-1], axis=1),axis=1)
            else:
                print("should never get here....")

        # Convert density profiles to N2
        for key in self.rr:
            N2 = (-g/case.rho0)*(np.diff(self.rr[key])/np.diff(self.zc)) # N2 at interior f points
            self.N2[key] = np.hstack((N2[0], N2, N2[-1]))                # Copy ends for bot/surf f points
            if np.any(N2<0):
                print("Unstable N2 profile detected")


class Projection(object):
    """
    Class to compute and store a projection. Handles u,v,w and p.
    """
    def __init__(self, VM, u=None, v=None, w=None, p=None, residual=False, correct=False, icut=None):
        """
        Constructor: Perform and store the projection
                     If 'residual' is True, we compute the residual
                     If 'correct' is True, we correct the w projection
                     by solving a linear system after the simple projection
                     If 'icut' is specified, we project only between y indices
                     icut and write zeros elsewhere (accelerates projection)
        """
        self.residual, self.correct, self.icut = residual, correct, icut
        self.um, self.ur = self.projwrapper(VM, u)
        self.vm, self.vr = self.projwrapper(VM, v)
        self.wm, self.wr = self.projwrapper(VM, w, wflag=True)
        self.pm, self.pr = self.projwrapper(VM, p)

    def projwrapper(self,VM,vv,wflag=False):
        """
        Wrapper function to call projection and residual
        """
        if vv is None:
            return None, None
        vvr=None
        if self.icut is None:
            vvm = self.project(VM, vv, wflag)
            if self.residual:
                vvr = self.computeresidual(VM, vv, vvm, wflag)
        else:
            i1, i2 = np.max([0, self.icut[0]]), np.min([self.icut[1],vv.shape[1]])
            vvm = np.zeros([VM.nmodes+1, vv.shape[1], vv.shape[2]])
            vvm[:,i1:i2,:] = self.project(VM, vv[:,i1:i2,:], wflag)
            if self.residual:
                vvr = np.zeros(vv.shape)
                vvr[:,i1:i2,:] = self.computeresidual(VM, vv[:,i1:i2,:], vvm[:,i1:i2,:], wflag)
        return vvm, vvr

    def project(self, VM, data, wflag):
        """
        Project 3D array data onto vertical modes. Returns a 3D array
        indexed as [m, y, x], where m=0 is the barotropic mode, m=1,2,3,... are
        baroclinic modes. If self.residual is True, we also return the residual.
        Lastly, wflag specifies to project onto dphic/N2c rather than phic,
        and self.correct controls applying the correction to w projection.
        """
        def projwrap(data, f, dz):  # Helper function to wrap this common code
            fdz = np.einsum('zm,z->zm',f,dz)
            try:
                import lpolib.fast_project as fp  # OpenMP accelerated C based projector
                return fp.project(data,fdz)
            except:
                return np.einsum('zyx,zm->myx',data,fdz)  # Default method

        dz = np.diff(VM.zf)/VM.H
        if wflag: # data is w
            # Estimate by ignoring cross terms
            pdata0 = projwrap(data,VM.dphic,dz)
            if self.correct:
                # Apply correction for nonzero cross terms
                A=np.einsum('zm,zn,z->mn',VM.dphic, VM.dphic, dz/VM.N2c)
                return np.einsum('nm,myx->nyx',np.linalg.inv(A),pdata0)
            else:  # No correction, just scale by c
                return np.einsum('myx,m->myx',pdata0,VM.c**2)
        else:     # data is u, v or p
            return projwrap(data,VM.phic,dz)

    def computeresidual(self, VM, data, pdata, wflag):
        """
        Computes the residual. The reconstruction is equal to (data - resid)
        """
        if wflag: # we're reconstructing w
            resid = data - np.einsum('myx,zm,z->zyx',pdata,VM.dphic,1/VM.N2c)
        else:     # u, v or p
            resid = data - np.einsum('myx,zm->zyx',pdata,VM.phic)
        return resid

class EnergyBreakdown(object):
    """
    Class to compute and hold the energy breakdown
    """
    def __init__(self, case=None, VM=None, P=None, u=None, v=None, cname=None, tmask=None):
        """
        Constructor
        - If only given cname, we load results from cname.
        - If given VM,P,u,v, we call computebreakdown(). If also
          given cname, we save results to cname
        """
        if all([x is None for x in [case,VM,P,u,v]]):
            # Load results from cname
            if os.path.isfile(cname):
                npzfile = np.load(cname)
                self.enn = npzfile['enn']
                self.emn = npzfile['emn']
                self.et,self.ep,self.er,self.vol = npzfile['etprv']
            else:
                print("Cache file " + cname + "does not exist, skipping...")
        elif all([x is not None for x in [case,VM,P,u,v]]):
            self.computebreakdown(case,VM,P,u,v,tmask)   # do the computation, save result
            if cname is not None:                  # save result
                np.savez(cname, enn=self.enn, emn=self.emn, etprv=[self.et,self.ep,self.er,self.vol])

    def computebreakdown(self,case,VM,P,u,v,tmask):
        """
        Computes the horizontal kinetic energy breakdown.
        Total energy is et = ep + er, and ep = sum(en) + sum(enm).
        """
        xT,yT = case.getXY('T');
        zf    = case.getZ('w')        # Nominal z grid at w points, 1D array
        dA    = (xT[0,1]-xT[0,0])*(yT[1,0]-yT[0,0]) # assume constant grid spacing
        dz    = np.diff(zf)
        self.emn = np.zeros([VM.nmodes+1,VM.nmodes+1])
        if tmask is None:
            self.tmask=1                 # default mask for entire domain
        else:
            self.tmask=tmask[:,1:-1,:]   # mask for integrating regions
        for mi in range(VM.nmodes+1):
            phim = VM.getmode(mi)[1]
            for ni in range(VM.nmodes+1):
                phin = VM.getmode(ni)[1]
                umun = P.um[mi,:,:]*P.um[ni,:,:]
                vmvn = P.vm[mi,:,:]*P.vm[ni,:,:]
                umunr = umun[1:-1,0:-1] + umun[1:-1,1:]   # 2umun at rho points
                vmvnr = vmvn[0:-1,1:-1] + vmvn[1:,1:-1]   # 2vmvn at rho points
                self.emn[mi,ni] = 0.25*np.sum(self.tmask*(umunr+vmvnr))*dA*np.sum(phin*phim*dz)
        self.enn = np.diag(self.emn)             # modal energies (vector)
        self.emn = self.emn - np.diag(self.enn)  # cross terms    (matrix, zeros on main diagonal)
        # For u = up + ur, energy is 0.5*u^2 = 0.5*up^2 + upur + 0.5*ur^2
        self.et = 0.5*volint(u**2       , v**2       , dA, dz, self.tmask)  # 0.5*(u^2 + v^2)
        self.ep = 0.5*volint((u-P.ur)**2, (v-P.vr)**2, dA, dz, self.tmask)  # 0.5*(up^2 + vp^2)
        self.er = self.et - self.ep                             # (upur + vpvr) + 0.5*(ur^2 + vr^2)
        self.vol = np.sum(self.tmask)*dA*case.H                 # volume


class EnergyHolder(object):
    def __init__(self, rname, vlist, tlist, rlist):
        self.rname = rname
        self.vlist, self.tlist, self.rlist = vlist, tlist, rlist
        self.nv,    self.nt,    self.nr    = len(vlist), len(tlist), len(rlist)
        self.E  = ndlist([self.nv,self.nt,self.nr])
    def name2index(self, vname, tname, rname):
        vi=[xi for xi,x in enumerate(self.vlist) if self.vlist[xi] == vname]
        ti=[xi for xi,x in enumerate(self.tlist) if self.tlist[xi] == tname]
        ri=[xi for xi,x in enumerate(self.rlist) if self.rlist[xi] == rname] # region name, not run name
        return vi[0],ti[0],ri[0]
    def set(self, E, vname, tname, rname):
        vi, ti, ri = self.name2index(vname, tname, rname)
        self.E[vi][ti][ri]=E
    def get(self, vname, tname, rname):
        vi, ti, ri = self.name2index(vname, tname, rname)
        return self.E[vi][ti][ri]
    def getts(self, vname, rname, item):
        vi, ti, ri = self.name2index(vname, self.tlist[0], rname)
        if item=="et":
            d= [self.E[vi][x][ri].et for x in range(self.nt)]
        if item=="ep":
            d= [self.E[vi][x][ri].ep for x in range(self.nt)]
        if item=="er":
            d= [self.E[vi][x][ri].er for x in range(self.nt)]
        if isinstance(item,int): # assume we specified a mode index here
            d= [self.E[vi][x][ri].enn[item] for x in range(self.nt)]
        return self.tlist, np.array(d)

def doprojection(N2p,reg,nmodes,u,v,w=None,p=None,residual=False):
    """
    Projection wrapper
     - N2p is a N2profile class containing N2 profiles
     - regs is a region key apready set in N2p (see N2profile class)
     - nmodes is number of modes to project
     - u,v are velocity data from loaduv
     - p is optional pressure data from loadp
     - residual boolean controls computation of residual
    """
    VM = VModes(N2p.zc,N2p.zf,N2p.N2[reg],nmodes)     # VM contains the eigenvalue problem solution
    P  = Projection(VM,u,v,w=w,p=p,residual=residual) # P contains u,v,[w],[p] and their projection/residual
    return VM,P


def yregmask(case,yreg,v):
    """
    Converts the yregions into masks with dimensions (1,ny,1) for broadcasting use
    """
    nreg = len(yreg)
    mask = [None]*nreg
    yT=case.getXY('T')[1]; yT=yT[:,0]
    yV=case.getXY('V')[1]; yV=yV[:,0]
    pt = case.hgrid.get_htype(v)
    if nreg == 1:
        if pt=="RHO" or pt=="U":
            mask[0]=np.ones_like(yT)[np.newaxis,:,np.newaxis]
        elif pt=="V" or pt=="PSI":
            mask[0]=np.ones_like(yV)[np.newaxis,:,np.newaxis]
        return mask
    for ri in range(nreg):
        ymin,ymax = yreg[ri]
        if pt=="RHO" or pt=="U":
            ii = np.where( (yT >  ymin) & (yT < ymax))[0] # cell centres in y
            if ii[0]==1:
                ii = np.hstack([ii[0]-1, ii])  # Include south boundary
            elif ii[-1] == len(yT)-2:
                ii = np.hstack([ii, ii[-1]+1]) # Include north boundary
            mask[ri]=np.zeros_like(yT)[np.newaxis,:,np.newaxis]
        elif pt=="V" or pt=="PSI":
            ii = np.where( (yV >= ymin) & (yV <= ymax))[0]   # indices for v
            mask[ri]=np.zeros_like(yV)[np.newaxis,:,np.newaxis]
        mask[ri][:,ii,:] = 1
    return mask


def mergeregions(case,yreg,v):
    """
    Merges data from yregions into one field via masking & summing
    """
    if isinstance(v,tuple): # return a tuple if f is a tuple
        return tuple([mergeregions(case,yreg,x) for x in v])
    if v is None:
        return v
    if len(v)==1:
        return v[0]
    mask = yregmask(case,yreg,v[0])
    nreg = len(yreg)
    vo = np.zeros(v[0].shape)
    for mi in range(nreg-1):
        ie=np.nonzero(mask[mi]+mask[mi+1]==2)[1]
        if np.any(ie):                 # average points that are
            mask[mi][:,ie,:]=0.5       # found on the boundary of
            mask[mi+1][:,ie,:]=0.5     # two adjacent regions
    for ri in range(nreg):
        vo += v[ri]*mask[ri]
    return vo


def eb_cachename(savepath,rname,vtype,ti,ri=None):
    # Helper function to build the complete energy breakdown cache filename
    cpath=os.path.join(savepath, rname)
    if not os.path.isdir(cpath):
        os.makedirs(cpath)
    if ri is None:
        fn = "EB.vtype." + vtype + ".ti." + str(ti) + ".npz"
    else:
        fn = "EB.vtype." + vtype + ".ti." + str(ti) + ".reg." + str(ri) + ".npz"
    return os.path.join(cpath, fn)


def avgden_cachename(savepath, rname, tilist):
    # Helper function to build the complete cache filename
    cpath=os.path.join(savepath, rname)
    if not os.path.isdir(cpath):
        os.makedirs(cpath)
    fn = "AvgDen.ti." + str(tilist[0]) + "-" + str(tilist[-1]) + ".npz"
    return os.path.join(cpath, fn)

def proj_cachename(savepath, rname, ti, reg, vtype):
    # Helper function to build the complete cache filename
    cpath=os.path.join(savepath, rname)
    if not os.path.isdir(cpath):
        os.makedirs(cpath)
    if type(reg) is str:
        regstr=reg
    if type(reg) is tuple:
        regstr=str(reg[0])+"-"+str(reg[1])

    fn = "Proj.ti." + str(ti) + ".reg." + regstr + ".type." + vtype + ".npz"
    return os.path.join(cpath, fn)

def ndlist(shape):
    """
    ndlist: generates an empty n-dimensional list, shape specifies dimensions
    """
    if shape:
        return [ndlist(shape[1:]) for i in range(shape[0])]

def geticut(case,creg):
    """
    Given a region, find y indices slightly below and slightly above the boundaries
    """
    yV=case.getXY('V')[1][:,0]
    jlow=np.max([0,np.where(yV>creg[0])[0][0] -4])
    jhigh=np.min([yV.shape[0], np.where(yV<=creg[1])[0][-1] +4])
    return jlow,jhigh

def skill(case,reg,ri,ve,v):
    """
    Skill metric used in the paper.
    Normalized mean square error over centre region subtracted from unity
    """
    mask=yregmask(case,reg,v)[ri].astype('bool')[0,:,0]
    return 1 - np.mean( (v[mask,:] - ve[mask,:])**2) / np.mean(ve[mask,:]**2)
