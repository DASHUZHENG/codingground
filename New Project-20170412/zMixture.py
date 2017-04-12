# -*- coding: utf-8 -*-
version="0.1"
import numpy
import sympy
import math

import scipy.integrate as integrate
import scipy.misc as diff

import zDatabase
import __future__   #the division calculation prolbem!

import logging  

from zChemical import zProperty
from zEquation import zEquation
from zEquation import zEquationMX

class zConcentration():
    def __init__(self,zmixture_example):
        self._concentration={}
        for content in zmixture_example:
            pass
        

class zMixture():
    
    def __init__(self,*args):
        
        self.logger=logging.getLogger("test.txt")
        
        self.vector=[]
        self.matrix=[]
        
        self.composition={}
        
        self.mole_content={}
        self.mass_content={}
        
        self.mass_conc={}
        self.mole_conc={}
        
        for arg in args:
            if isinstance(arg,zProperty):
                self.composition[arg.name]=arg
        
        self.vector=[x for x in self.composition]
        self.matrix=[(x,y) for x in self.vector for y in self.vector]
        
        
        
        self._eoslist=["IDEAL","RK","RKS","PR"]
        self._activitylist=["WILSON","NTRL","UNIQUAC"]
        self._defaultactivity={\
        "WILSON":(0,0,0,0,0,0,0,0,0,1000,0,0),\
        "NRTL":(0,0,0,0,0,0,0,0,0,0,0,1000),\
        "UNIQUAC":(0,0,0,0,0,0,0,0,0,0,0,1000)}
        
        
        
        self._thermoset=['PHIVMX','PHILMX','HVMX','HLMX','GVMX','GLMX','SVMX','SLMX','VVMX','VLMX','PHISMX','HSMX','GSMX','SSMX','VSMX','WSLMX','HCSLMX']
        self._thermosubset=['PHILPCMX','DHVMX','DHLMX','DGVMX','DGLMX','DSVMX','DSLMX','PHISOCMX']
        
        __slot__=[]
        
        #self.conc_setting()
        
    def eos_set(self,eos="RK",activity="WILSON",database="zPro.db",location="biWILSON"):
        
        if eos in self._eoslist:
            self.eos=eos
            #!Database Problem
            self.eosset=(eos,self.tc,self.pc,self.vc,self.omega)
            #Is this enough for EOS tuple?
        else:
            self.eos="RK"
            self.eosset=(eos,self.tc,self.pc,self.vc,self.omega)
            #Is this enough for EOS tuple?
            print "Default Equation is RK"
        
        if activity in self._activitylist:
            
            self.activity=activity
            self.gamma_method={}
            
            calcfilter=lambda x: isinstance(x,float)
            
            for comp in self.matrix:
            
                try:
                    pair1=self.matrix[comp]
                    pair2=(pair1[1],pair1[0])
                    
                    search_activity=db.Locate(database,location,pair1,"*","Field1")
                    
                    if search_activity==None:
                        search_activity=db.Locate(database,location,pair2,"*","Field1")
                    
                    if search_activity==None:
                        search_activity=self._defaultactivity[self.activity]
                    
                    #!!All the prime key change to "Field1"!!!
                    #print("Debug tdep property load member %s %s" % (btt,result2))
                    
                    filtered_activity=filter(calcfilter,search_activity)
                    activity_coef=list(filtered_activity)
                    activity_coef.insert(0,self.activity)
                    self.gamma_method[comp]=tuple(activity_coef)
                
                except Exception,error:
                    self.logger.info("Activilty Load Error")
                    #print("KeyWord:%s;\nTdep_property_database error: %s\n" % (att,error))
                    self.gamma_method[comp]=None
                    
                    #!Not Finish Yet
    
   
   
   
    def eos_setting(self):
        pass
    
    
    def Gamma(self,temp,pres,t0=298.15,p0=100000,route="1"):
        gamma_list={}
        if route=="1":
            for comp in self.vector:
                gamma_list[comp]=1
                
        elif route=="WILSON":
            pass
            
        pass
        
        return gamma_list
        
    def Henry(self,temp,pres,t0=298.15,p0=100000,route="Henry"):
        pass
    
    
    
    def Vvmx(self,temp,pres,t0=298.15,p0=100000,route="RK"):
        
        if route=="1":
            pass
        
        elif route=="2":
            pass
        
        elif route=="RK":
            
            equation=zEquationMX()
            
            try:
                result=getattr(equation,"RKMX")(self.vector,self.matrix,self.composition,self.mole_conc,temp,pres,0,tag="V")
                
                #print("Vv Debug %s %s %s" %(self.eosset,result,8.3145*temp/pres))
            except Exception, error:
                print ("KeyWord:%s + %s\nMistake:%s\n"%("Vvmx","Test",error))
                result= 0
        else:
            
            pass
            
            
        return result
        
    def Vlmx(self,temp,pres,t0=298.15,p0=100000,route="2"):
        
        if route=="1":
            pass
        elif route=="2":
            vl_list={}
            vlmx=0
            
            for comp in self.vector:
                vl_list[comp]=self.composition[comp].Vl(temp,pres)
                vlmx=vlmx+vl_list[comp]*self.mole_conc[comp]
            
        result=vlmx
        
        return result
            
    def Phivmx(self,temp,pres,t0=298.15,p0=100000,route="RK"):
        
        if route=="1":
            pass
        
        elif route=="RK":
            
            equation=zEquationMX()
            
            try:
                result=getattr(equation,"RKMX")(self.vector,self.matrix,self.composition,self.mole_conc,temp,pres,0,tag="PHIV")
                
                #print("Vv Debug %s %s %s" %(self.eosset,result,8.3145*temp/pres))
            
            except Exception, error:
                print ("KeyWord:%s + %s\nMistake:%s\n"%("Phivmx","Test",error))
                result= 0
                #!Not Yet Finished
        
        return result
    
    def Philmx(self,temp,pres,t0=298.15,p0=100000,route="2"):
        
        if route=="1":
            pass
        
        if route=="2":
            #Default Lewis-Randall or Rault Law
            
            gamma_list=self.Gamma(temp,pres,t0,p0,route="1")
            
            for comp in self.vector:
                phil_list[comp]=self.zcomposition[comp].Phil(temp,pres)
                philmx_list[comp]=phil_list[comp]*gamma_list[comp]
        
        
        return philmx_list
                
                
            
            
        
        
        

        
if __name__=="__main__":
    #a initiated
    a=zProperty("H2O")
    #b initiated
    b=zProperty("H2O")
    b.name="ACETIC"
    
    test=zMixture(a,b)
    print test.vector
    print test.matrix
    print test.composition

    test.mole_content["H2O"]=5
    test.mole_content["ACETIC"]=5
    test.mole_conc["H2O"]=0.5
    test.mole_conc["ACETIC"]=0.5
    
    print test.mole_content
    
    print "Vvmx",test.Vvmx(500,300000)
    print "Vlmx",test.Vlmx(298.15,100000)
    print "Phiv",test.Phivmx(500,30000)
    print "Phil",test.Phivmx(298.15,100000)

