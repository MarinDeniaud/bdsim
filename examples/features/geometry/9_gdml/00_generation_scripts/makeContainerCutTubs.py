import os as _os
import pyg4ometry.gdml as _gd
import pyg4ometry.geant4 as _g4
import pyg4ometry.visualisation as _vi

normal     = 1 
flat_ends  = 2

def Test(vis = False, interactive=False, type = normal,n_slice=16) :
    reg = _g4.Registry()
    
    # defines 
    wx = _gd.Constant("bx","1",reg,True)
    wy = _gd.Constant("by","1",reg,True)
    wz = _gd.Constant("bz","1",reg,True)

    # pi         = _gd.Constant("pi","3.1415926",reg,True)
    ctrmin     = _gd.Constant("trmin","0",reg,True)
    ctrmax     = _gd.Constant("trmax","10.0",reg,True)
    ctz        = _gd.Constant("tz","50",reg,True)
    ctstartphi = _gd.Constant("startphi","0",reg,True)
    ctdeltaphi = _gd.Constant("deltaphi","2.0*pi",reg,True)
    ctlowx     = _gd.Constant("ctlowx","-1",reg,True)
    ctlowy     = _gd.Constant("ctlowy","-1",reg,True)
    ctlowz     = _gd.Constant("ctlowz","-1",reg,True)
    cthighx    = _gd.Constant("cthighx","1",reg,True)
    cthighy    = _gd.Constant("cthighy","1",reg,True)
    cthighz    = _gd.Constant("cthighz","1",reg,True)
            
    wm = _g4.Material(name="G4_Galactic") 
    bm = _g4.Material(name="G4_Fe") 

    # solids
    bs = _g4.solid.Box("ws",wx,wy,wz, reg,"mm")
    ws = _g4.solid.CutTubs("ts",ctrmin,ctrmax,ctz,ctstartphi,ctdeltaphi,[ctlowx,ctlowy,ctlowz],[cthighx,cthighy,cthighz],reg,"mm","rad",nslice=n_slice)
        
    # structure 
    wl = _g4.LogicalVolume(ws, wm, "wl", reg)
    ctl = _g4.LogicalVolume(bs, bm, "ctl", reg)
    ctp = _g4.PhysicalVolume([0,0,0],[0,0,0],  ctl, "ct_pv1", wl, reg) 
    
    # set world volume
    reg.setWorld(wl.name)
    
    # gdml output 
    w = _gd.Writer()
    w.addDetector(reg)
    w.write(_os.path.join(_os.path.dirname(__file__), "../gdmls/container-cuttubs.gdml"))
    w.writeGmadTester(_os.path.join(_os.path.dirname(__file__),"../09_container_cuttubs.gmad"),"gdmls/container-cuttubs.gdml")

    # visualisation
    v = None
    if vis : 
        v = _vi.VtkViewer()
        v.addLogicalVolume(reg.getWorldVolume())
        v.view(interactive=interactive)

if __name__ == "__main__":
    Test()
