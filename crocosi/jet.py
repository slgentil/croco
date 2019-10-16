
# This piece of code needs to be reworked

def getini(self, config=1):
    """
    getini: return jet initial conditions
    """

    # init grid, grid will be updated later on
    u=np.zeros((self.N,self.M+1,self.L))
    v=np.zeros((self.N,self.M,self.L+1))
    zeta=np.zeros((self.M+1,self.L+1))
    z_r=self.getZ('RHO', zeta=zeta)
    z_w=self.getZ('W', zeta=zeta)

    # a building function
    def asymth(zz,zasym,dzasym):
        return (zz-zasym)*(1.0+0.5*((zz-zasym)+np.abs(zz-zasym))/dzasym) \
                +zasym

    # from ana_initial.F
    # southern profile
    rhomaxs=27.75;
    bgdrhos=9.8e-6;
    zs1=-1000;
    dzs=700;
    # drhos=1.4;
    # northern profile
    rhomaxn=27.7573;
    bgdrhon=9.8e-6;
    zn1=-400;
    dzn=300;
    # surface Charney mode
    # drhosfs=1.5;
    # drhosfn=0.0;
    # z0=-300;
    drhosf=0.00;
    # z0p=-110;
    # aponte jet
    # z0p=z0;
    # alpha1=0.0075;

    if config==1:
        # cfg1 ~ setup 1 (albeit for alpha1 and z0p)
        # strong charney mode
        print('Configuration 1: strong charney mode (~ setup 1)')
        drhos=1.4;
        drhosfs=1.5;
        drhosfn=0.0;
        alpha1=0.0075;
        alpha2=1.;
        z0=-300;
        z0p=z0;
        #
        drhosfs1=drhosfs*alpha1;
        drhosfs2=drhosfs*alpha2;
        drhosfn1=drhosfn*alpha1;
        drhosfn2=drhosfn*alpha2;
    elif config==2:
        # cfg2 ~ setup 2 (albeit for alpha1
        # Phillips, weak charney mode
        print('Configuration 2: Phillips, weak charney mode (~ setup 2)')
        drhos=1.4;
        drhosfs=1.5;
        drhosfn=0.0;
        alpha1=0.0075;
        alpha2=0;
        z0=-300;
        z0p=z0;
        #
        drhosfs1=drhosfs*alpha1;
        drhosfs2=drhosfs*alpha2;
        drhosfn1=drhosfn*alpha1;
        drhosfn2=drhosfn*alpha2;
    elif config==3:
        # cfg 3
        # setup 1: Phillips>Charney
        print('Configuration 3: Phillips>Charney (setup 1)')
        drhos=1.4;
        drhosfs1=0; # drhosf
        drhosfn1=0; # drhosf
        drhosfs2=1.5; # drhosfs
        drhosfn2=0.0; # drhosfn
        #drhosfs=0.0; # test
        #drhosfn=1.5; # test
        z0=-300;
        z0p=-110;
    elif config==4:
        # cfg 4
        # setup 2: Phillips
        print('Configuration 4: Phillips (setup 2)')
        drhos=1.4;
        drhosfs1=0.; # drhosf
        drhosfn1=0.; # drhosf
        drhosfs2=0.; # drhosfs
        drhosfn2=0.; # drhosfn
        z0=-300;
        z0p=-110;
    elif config==5:
        # cfg 5
        # setup 4: Phillips + surface temperature anomaly
        print('Configuration 5: Phillips + surface temperature anomaly (setup 4)')
        drhos=1.4;
        drhosfs1=0.71; # drhosf
        drhosfn1=0.71; # drhosf
        drhosfs2=1.5; # drhosfs
        drhosfn2=0.; # drhosfn
        z0=-300;
        z0p=-110;
    elif config==6:
        # cfg 6
        # setup 5: strong charney
        print('Configuration 6: strong charney (setup 5)')
        drhos=0.;
        drhosfs1=0.; # drhosf
        drhosfn1=0.; # drhosf
        drhosfs2=2.5; # drhosfs
        drhosfn2=1.0; # drhosfn
        z0=-600;
        z0p=-110;

    # zonal perturbation
    cff_perturb=0.02;

    ### ! --- Build northern and southern density profiles ---

    # ! First compute background density profiles associated with
    # ! a gentle stratification that does not depend
    # ! on jet side. It is there to ensure static stability.

    i=0;j=0; # zk = z_r(istr,jstrR,k)
    h = self.H

    rhoprof = np.zeros((self.N,2))
    rhoprof[:,0]=rhomaxs-bgdrhos*(z_r[:,j,i]+h);
    rhoprof[:,1]=rhomaxn-bgdrhon*(z_r[:,j,i]+h);
    rho0=rhoprof[:]

    # ! Second, get main north/south contrast with a distorded
    # ! tanh shape of variable amplitude.Distorsion of the tanh
    # ! is done with the asym function that increases
    # ! stratification in the upper ocean

    dzs_a=1.3*dzs;
    zk=asymth(z_r[:,j,i],zs1,dzs_a);
    rhoprof[:,0]= rhoprof[:,0]-drhos*(0.5+0.5*np.tanh((zk-zs1)/dzs));

    dzn_a=1.3*dzn;
    drhon=-(rhoprof[self.N-1,0]-rhoprof[self.N-1,1])/(0.5+0.5*np.tanh((zk[-1]-zn1)/dzn));
    zk=asymth(z_r[:,j,i],zn1,dzn_a);
    rhoprof[:,1]= rhoprof[:,1]-drhon*(0.5+0.5*np.tanh((zk-zn1)/dzn));

    rho1 = rhoprof[:];

    zk=z_r[:,j,i];
    rhoprof[:,0] += -drhosfs1*(np.exp((zk-z0p)/np.abs(z0p)))/(np.exp(1.));
    #      -drhosfs*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.));
    rhoprof[:,0] += -drhosfs2*0.5*(1+np.tanh((zk-z0)/np.abs(z0)))/np.tanh(1.);
    #      -drhosfs*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.);
    rhoprof[:,1] += -drhosfn1*(np.exp((zk-z0p)/np.abs(z0p)))/(np.exp(1.));
    #      -drhosfn*alpha1*(exp((zk-z0p)/abs(z0p)))/(exp(1.));
    rhoprof[:,1] += -drhosfn2*0.5*(1+np.tanh((zk-z0)/np.abs(z0)))/np.tanh(1.);
    #      -drhosfn*alpha2*0.5*(1+tanh((zk-z0)/abs(z0)))/tanh(1.);

    # averages the profiles in order to have smaller difference
    # and ultimately a weaker jet
    print("Jet ywidth is ", self.jet_ywidth/1e3, "km")
    print("Jet weight is ", self.jet_weight)
    rhoprof[:,1] = self.jet_weight * rhoprof[:,1] \
                    + ( 1 - self.jet_weight) * rhoprof[:,0];

    rho2 = rhoprof[:];

    # horizontal indices
    x = (np.arange(-1,self.Lm+1)+0.5)/(self.Lm-2.0)
    y = (np.arange(-1,self.Mm+1)+0.5)/self.Mm -0.5
    y = np.tile(y.reshape((1,y.size,1)),(self.N,1,self.L+1))
    x = np.tile(x.reshape((1,1,x.size)),(self.N,self.M+1,1))


    # if flag_jet_perturb:
    y = y + cff_perturb*np.exp(z_r/1000.) \
                      *np.exp(-(x-0.5)**2/0.05) \
                      *( 0.5*np.sin(2.*np.pi*x) \
                        +0.5*np.sin(6.*np.pi*x) );

    y = y*np.pi*self.hgrid.el/self.jet_ywidth + np.pi/2;

    Fyz=0.5-(y-np.sin(y)*np.cos(y))/np.pi;
    Fyz[np.where(y<0)]=0.5;
    Fyz[np.where(y>np.pi)]=-0.5;
    Fyz[:]=Fyz[:]+0.5;

    rhos = np.tile(rhoprof[:,0].reshape((self.N,1,1)),(1,self.M+1,self.L+1))
    rhon = np.tile(rhoprof[:,1].reshape((self.N,1,1)),(1,self.M+1,self.L+1))
    rho = Fyz*rhos + (1.-Fyz)*rhon;

    # compute pressure and pressure gradients
    flag_flux=False;
    dpdx,dpdy,P = get_gradp(rho,z_r,z_w,self,flag_flux);

    # adjust sea level in order to have zero pressure signal at the bottom
    zeta = -P[0,:,:]/9.81; # P is already divided by rho0
    zeta += -zeta[:].mean();

    # recompute depth and pressure gradients
    z_r=self.getZ('RHO', zeta=zeta)
    z_w=self.getZ('W', zeta=zeta)
    dpdx,dpdy,P = get_gradp(rho,z_r,z_w,self,flag_flux);

    # get Coriolis parameter
    f=self.hgrid.f[:]

    # compute geostrophic velocities
    #print dpdx.shape, dpdy.shape
    u[:,1:-1,:] = (dpdy[:,:-1,:-1]+dpdy[:,:-1,1:]+ \
                   dpdy[:,1:,:-1]+dpdy[:,1:,1:])*.25;
    u[:,0,:] = u[:,1,:]
    u[:,-1,:] = u[:,-2,:]
    u = u/f[:,:-1]
    #u = np.concatenate((u[:,[0],:],u,u[:,[-1],:]),axis=1)/f[:,:-1];
    v[:,:,1:-1] = (dpdx[:,:-1,:-1]+dpdx[:,:-1,1:]+ \
         dpdx[:,1:,:-1]+dpdx[:,1:,1:])*.25;
    # assumes zonal periodicity
    v[:,:,0]=v[:,:,-2]
    v[:,:,-1]=v[:,:,1]
    #print v.shape, f.shape
    #v = np.concatenate((v[:,[0],:],v),axis=1)/f[:,:-1];

    return zeta, rho, u, v
