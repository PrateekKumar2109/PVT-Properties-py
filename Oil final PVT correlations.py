def Rs_gas_solubility(o_api,temp,Pb_pressure,sg,press_sep,temp_sep):
  #First correlation is Standings (1981) based on California samples#
  #Second correlation is Vasquez-Beggs (1980) #
   #Third correlation is Glaso (1980) based on North Sea samples#
 #Fourth correlation is Marhoun (1988) based on Middle East samples#
 #Fifth correlation is Petrosky-Farshad (1993) based on Gulf of Mexico samples#   
    import numpy as np
    import math
    temp=(np.array(temp)*1.8)+32 #Temp in deg C conversion to deg F#
    temp=temp+459.67
    temp_sep=np.array(temp_sep)+459.67
    pressure=np.array(Pb_pressure)+14.23 #Pb pressure in psig#
    press_sep=np.array(press_sep)+14.23
    #First Standing Correlation#
    x=(0.0125*np.array(o_api))-(0.00091*(np.array(temp)-460))#temp in rankine#
    Rs_standing=(np.array(sg))*((((np.array(pressure)/18.2)+1.4)*10**x)**1.2048)
    Rs_standing=Rs_standing*0.1781 
    #sg gas gravity at seperator pressure & tmperature#
    
    #Second Vasquez-Beggs Correlation#
    sg_ref_sep=[]
    Rs_beggs=[]
    for i in range(len(o_api)):
     if o_api[i]>30:
        c1=0.0178;c2=1.1870;c3=23.931
     else:
        c1=0.0362;c2=1.0937;c3=25.7240
     sg_ref_sep1=(sg[i])*(1+(0.00005912*(o_api[i])*((temp_sep[i])-460)*np.log10((press_sep[i])/114.7)))
     rs_beggs=c1*sg_ref_sep1*(pressure[i]**c2)*np.exp(c3*(o_api[i])/(temp[i]))
     rs_beggs=rs_beggs*0.1781   #Conversion from scf/STB to v/v#
     sg_ref_sep.append(sg_ref_sep1)
     Rs_beggs.append(rs_beggs) 
        
    #Third Glaso Correlation#
     pb_x=np.array(10**((2.8869)-(14.1811-3.3093*(np.log10(np.array(pressure))))**0.5))
     Rs_glaso=(np.array(sg))*((pb_x*((np.array(o_api)**0.989)/((np.array(temp)-460)**0.172)))**1.2255)
     Rs_glaso=Rs_glaso*0.1781  
        
    #Fourth Marhouns Correlation#
     a=185.843208;b=1.877840;c=-3.1437;d=-1.32657;e=1.398441
     Rs_marhoun= (a*(np.array(pressure))*(np.array(sg)**b)*((141.5/(131.5 +np.array(o_api)))**c)*(np.array(temp)**d))**e
     Rs_marhoun=Rs_marhoun*0.1781
    #fifth Petrosky-Farshad Correlation
     x_1=(0.0007916*(np.array(o_api)**1.5410))-(0.00004561*((np.array(temp)-460)**1.3911))
     Rs_petrosky=((10**x_1)*(np.array(sg)**0.8439)*(12.340+(np.array(pressure)/112.727)))**1.73184
     Rs_petrosky=Rs_petrosky*0.1781
        
    return( Rs_standing,Rs_beggs,Rs_glaso,Rs_marhoun,Rs_petrosky)


def Bubble_point_pressure(o_api,temp,Pb_pressure,sg,press_sep,temp_sep,Rs):
  #First correlation is Standings (1981) based on California samples#
  #Second correlation is Vasquez-Beggs (1980) #
   #Third correlation is Glaso (1980) based on North Sea samples#
 #Fourth correlation is Marhoun (1988) based on Middle East samples#
 #Fifth correlation is Petrosky-Farshad (1993) based on Gulf of Mexico samples#   
    import numpy as np
    import math
    Rs=np.array(Rs)/0.1781
    temp=(np.array(temp)*1.8)+32 #Temp in deg C conversion to deg F#
    temp=temp+459.67
    temp_sep=np.array(temp_sep)+459.67
    pressure=np.array(Pb_pressure)+14.23 #Pb pressure in psig#
    press_sep=np.array(press_sep)+14.23
    #First Standing Correlation#
    a=-(0.0125*np.array(o_api))+(0.00091*(np.array(temp)-460))#temp in rankine#
    Pb_standing=((((np.array(Rs)/np.array(sg))**0.83)*(10**a))-(1.4))*18.2
    #sg gas gravity at seperator pressure & tmperature#
    
    #Second Vasquez-Beggs Correlation#
    sg_ref_sep=[]
    Pb_beggs=[]
    for i in range(len(o_api)):
     if o_api[i]>30:
        c1=56.18;c2=0.84246;c3=10.393
     else:
        c1=27.624;c2=0.914328;c3=11.172
     a_1=-c3*o_api[i]/temp[i]
     sg_ref_sep1=(sg[i])*(1+(0.00005912*(o_api[i])*((temp_sep[i])-460)*np.log10((press_sep[i])/114.7)))
     pb_beggs=((c1*Rs[i]/sg_ref_sep1)*(10**a_1))**c2
     sg_ref_sep.append(sg_ref_sep1)
     Pb_beggs.append(pb_beggs) 
        
    #Third Glaso Correlation#
     a=0.816;b=0.172;c=-0.989
     pb_x=((np.array(Rs)/np.array(sg))**a)*((np.array(temp)-460)**b)*((np.array(o_api)**c))
     Pb_glaso= 10**(1.7669+(1.7447*np.log10(np.array(pb_x)))-(0.30218*((np.log10(np.array(pb_x)))**2)))
      
    #Fourth Marhouns Correlation#
     a=0.00538088;b=0.715082;c=-1.87784;d=3.1437;e=1.32657
     Pb_marhoun= (a*(np.array(Rs)**b)*(np.array(sg)**c)*((141.5/(131.5 +np.array(o_api)))**d)*(np.array(temp)**e))
    
    #fifth Petrosky-Farshad Correlation
     x_1=(0.0007916*(np.array(o_api)**1.5410))-(0.00004561*((np.array(temp)-460)**1.3911))
     Pb_petrosky=((112.727*(np.array(Rs)**0.577421))/((10**x_1)*(np.array(sg)**0.8439)))-1391.051
    
    return( Pb_standing,Pb_beggs,Pb_glaso,Pb_marhoun,Pb_petrosky)

#Oil FVF#
def oil_formation_volume_factor(o_api,temp,Pb_pressure,sg,press_sep,temp_sep,Rs,res_press):
  #First correlation is Standings (1981) based on California samples#
  #Second correlation is Vasquez-Beggs (1980) #
   #Third correlation is Glaso (1980) based on North Sea samples#
 #Fourth correlation is Marhoun (1988) based on Middle East samples#
 #Fifth correlation is Petrosky-Farshad (1993) based on Gulf of Mexico samples#   
    import numpy as np
    import math
    Rs=np.array(Rs)/0.1781
    temp=(np.array(temp)*1.8)+32 #Temp in deg C conversion to deg F#
    temp=temp+459.67
    temp_sep=np.array(temp_sep)+459.67
    pressure=np.array(Pb_pressure)+14.23 #Pb pressure in psig#
    press_sep=np.array(press_sep)+14.23
    res_press=np.array(res_press)+14.23
    Bob=1.176
    

    #First Standing Correlation#
    #temp in rankine#
    Bo_standing=((((((np.array(sg)/(141.5/(131.5 +np.array(o_api))))**0.5)*(np.array(Rs)))+((np.array(temp)-460)*1.25))**1.2)*0.000120)+0.9759
    #sg gas gravity at seperator pressure & temperature#
    
    #Second Vasquez-Beggs Correlation#
    sg_ref_sep=[]
    Bo_beggs=[]
    Bo_undersaturated_beggs=[]
    
    for i in range(len(o_api)):
     sg_ref_sep1=(sg[i])*(1+(0.00005912*(o_api[i])*((temp_sep[i])-460)*np.log10((press_sep[i])/114.7)))
     if res_press[i]<pressure[i]:
        if o_api[i]>30:
         c1=4.67*(10**-4);c2=1.1*(10**-5);c3=1.337*(10**-9)
         bo_beggs=(((o_api[i])/sg_ref_sep1)*((temp_sep[i])-520)*(c2+(c3*Rs[i])))+(c1*Rs[i])+1
         Bo_beggs.append(bo_beggs) 
        else:
         c1=4.677*(10**-4);c2=1.751*(10**-5);c3=-1.811*(10**-8)
         bo_beggs=(((o_api[i])/sg_ref_sep1)*((temp_sep[i])-520)*(c2+(c3*Rs[i])))+(c1*Rs[i])+1
         Bo_beggs.append(bo_beggs) 
      
     else: 
         a_b=(10**-5)*(-1433+(5*Rs[i])+(17.2*((temp[i])-460))-(1180*sg_ref_sep1)+(12.61*o_api[i]))
         bo_undersaturated_beggs=Bob*(np.exp(-a_b*(np.log(res_press[i]/pressure[i]))))
         Bo_undersaturated_beggs.append(bo_undersaturated_beggs)
     sg_ref_sep.append(sg_ref_sep1)
     
    #Third Glaso Correlation#
     b=((((np.array(sg)/(141.5/(131.5 +np.array(o_api))))**0.526)*(np.array(Rs)))+((np.array(temp)-460)**0.968))
     a=-(((np.log10(b))**2)*0.27683)+((np.log10(b))*2.91329)-6.58511
     Bo_glaso= (10**(a))+1.0
      
    #Fourth Marhouns Correlation#
     a=0.742390;b=0.323294;c=-1.202040
     f=(((np.array(sg))**b)*((141.5/(131.5 +np.array(o_api)))**c)*((np.array(Rs))**a))
     Bo_marhoun= (0.497069)+(np.array(temp)*(0.862963*(10**-3)))+(f*(0.182594*(10**-2)))+((f**2)*(0.318099*(10**-5)))
    
    #fifth Petrosky-Farshad Correlation
     x_1=(0.0007916*(np.array(o_api)**1.5410))-(0.00004561*((np.array(temp)-460)**1.3911))
     a_1=(((((np.array(sg))**0.2914)/((141.5/(131.5 +np.array(o_api)))**0.6265))*((np.array(Rs))**0.3738))+(((np.array(temp)-460)**0.5371)*0.24626))**3.0936
     Bo_petrosky=(7.2046*(10**-5)*a_1)+1.0113
     a_p=((np.array(sg))**0.1885)*((np.array(o_api))**0.6265)*((np.array(Rs))**0.69357)*((np.array(temp)-460)**0.6729)*(4.1646*(10**-7)) 
     Bo_undersaturated_petrosky=Bob*(np.exp(-a_p*((res_press**0.4094)-(pressure**0.4094))))
    return(Bo_standing,Bo_beggs,Bo_undersaturated_beggs,Bo_glaso,Bo_marhoun,Bo_petrosky,Bo_undersaturated_petrosky)

#Oil Viscosity#
def oil_viscosity(o_api,temp,Pb_pressure,sg,press_sep,temp_sep,Rs,res_press):
  #First correlation is Beals (1981) based on California samples#
  #Second correlation is Robinson (1980) #
   #Third correlation is Glaso (1980) based on North Sea samples#
 #Fourth correlation is Marhoun (1988) based on Middle East samples#
 #Fifth correlation is Petrosky-Farshad (1993) based on Gulf of Mexico samples#   
    import numpy as np
    import math
    Rs=np.array(Rs)/0.1781
    temp=(np.array(temp)*1.8)+32 #Temp in deg C conversion to deg F#
    temp=temp+459.67
    temp_sep=np.array(temp_sep)+459.67
    pressure=np.array(Pb_pressure)+14.7 #Pb pressure in psig#
    press_sep=np.array(press_sep)+14.7
    #First Standing Correlation#
    #temp in rankine#
    a=10**(0.43+(8.33/np.array(o_api)))
    Visc_dead_oil_standing=(((1.8*(10**7))/(np.array(o_api)**4.53))+0.32)*(360/(temp-260))
    e=3.74*(10**-3)*(np.array(Rs));d=1.1*(10**-3)*(np.array(Rs));c=8.62*(10**-5)*(np.array(Rs))
    b=(0.68/(10**c))+ (0.25/(10**d))+(0.062/(10**e))    
    a_1= (np.array(Rs))*((2.2*(10**-7)*(np.array(Rs)))-(7.4*(10**-4)))                                                                          
    Visc_saturated_oil_standing=(Visc_dead_oil_standing**b)*(10**a_1)
    #sg gas gravity at seperator pressure & temperature#
    
    #Second Robinson-Beggs Correlation#
    
    y=10**(3.0324-(0.02023*np.array(o_api)))
    x=y*((np.array(temp)-460)**-1.163)
    Visc_dead_oil_beggs=(10**x)-1
    a_2= ((np.array(Rs)+100)**-0.515)*10.715 
    b_2= ((np.array(Rs)+150)**-0.338)*5.44                                                                             
    Visc_saturated_oil_beggs=(Visc_dead_oil_beggs**b_2)*(a_2) 
    a_3= (np.array(pressure)*(-3.9*(10**-5)))-5
    m=2.6*(10**a_3)*(np.array(pressure)**1.187)
    Visc_undersaturated_oil_beggs =Visc_saturated_oil_beggs*(((np.array(res_press))/(np.array(pressure)))**m)                                                                            
    
    return(Visc_dead_oil_standing, Visc_saturated_oil_standing,Visc_dead_oil_beggs,Visc_saturated_oil_beggs,Visc_undersaturated_oil_beggs)
