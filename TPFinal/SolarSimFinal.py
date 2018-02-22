#Zach Pomper: zbp

import numpy as np
import _thread as thread
from numpy import float32,float64
import math, string, copy, random, datetime,os
from PIL import ImageTk,Image
from tkinter import *
    
def scaleVec(vec,data):
    x,y,z=vec[0],vec[1],vec[2]
    dist=x+data.windowDep
    scaler=data.windowDep*(1/dist) 
    maxz=data.windowDep*math.tan(data.vertViewTheta/180*math.pi)
    maxy=data.windowDep*math.tan(data.horizViewTheta/180*math.pi)
    return [x*scaler,y*scaler*data.width//2/maxy,z*scaler*data.height//2/maxz] 
    
def doMove(Forward,Right,Up,data): #I had already known how to do this, but
    #google 3-space matrix translation for explanations if curious
    moveinc=data.moveinc
    translator=np.array([moveinc*Forward,moveinc*Right,moveinc*Up])
    tranmat=np.array([(1,0,0,translator[0]),(0,1,0,translator[1]),
                      (0,0,1,translator[2]),(0,0,0,1)])
    temp=np.transpose(np.dot(tranmat,(np.vstack(
                     (np.transpose(data.vecList),np.ones(len(data.vecList)))))))
    homogs=temp[0:,3]
    homogs=np.transpose(np.vstack((homogs,homogs,homogs)))
    data.vecList=temp[0:,0:3]*homogs**(-1)
    movePlanetInfo(tranmat,data)

def movePlanetInfo(tranmat,data):
    for planet in Planet.planets:
        tempLastPos=planet.lastPos.reshape((3,1))
        tempMat=np.dot(tranmat,np.vstack(
                (tempLastPos,np.ones(1)))).reshape((1,4))[0]
        planet.lastPos=tempMat[0:3]*tempMat[3]**(-1)
        
def doRotate(roll,pitch,yaw,data): #pretty simple stuff
    toRad=math.pi/180
    psi=pitch*data.psiRotInc*toRad
    theta=roll*data.thetaRotInc*toRad
    gamma=yaw*data.gammaRotInc*toRad
    xangmat=np.array([
    (1,0,0),
    (0,math.cos(theta),-math.sin(theta)),
    (0,math.sin(theta),math.cos(theta))])
    yangmat=np.array([
    (math.cos(psi),0,math.sin(psi)),
    (0,1,0),
    (-math.sin(psi),0,math.cos(psi))])
    zangmat=np.array([
    (math.cos(gamma),-math.sin(gamma),0),
    (math.sin(gamma),math.cos(gamma),0),
    (0,0,1)])
    q=np.dot(np.dot(xangmat,yangmat),zangmat)
    data.vecList=np.dot(data.vecList,q)
    rotatePlanetInfo(roll,pitch,yaw,data)
    data.systemCenter=rotatevec(roll,pitch,yaw,data.systemCenter,data)
    normalizeCoords(data)
def rotatePlanetInfo(roll,pitch,yaw,data):
    for planet in Planet.planets:
        planet.vel=rotatevec(-roll,-pitch,-yaw,planet.vel,data) 
        planet.lastPos=rotatevec(-roll,-pitch,-yaw,planet.lastPos,data)

def rotatevec(roll,pitch,yaw,vec,data):
    toRad=math.pi/180
    psi=pitch*data.psiRotInc*toRad
    theta=roll*data.thetaRotInc*toRad
    gamma=yaw*data.gammaRotInc*toRad
    xangmat=np.array([
    (1,0,0),
    (0,math.cos(theta),-math.sin(theta)),
    (0,math.sin(theta),math.cos(theta))])
    yangmat=np.array([
    (math.cos(psi),0,math.sin(psi)),
    (0,1,0),
    (-math.sin(psi),0,math.cos(psi))])
    zangmat=np.array([
    (math.cos(gamma),-math.sin(gamma),0),
    (math.sin(gamma),math.cos(gamma),0),
    (0,0,1)])
    temp=np.dot(np.dot(xangmat,yangmat),zangmat)
    return np.dot(temp,vec)

def normalizeCoords(data):
    delta=data.systemCenter-data.systemCenterCopy
    data.moveBuffer=(
    [data.moveBuffer[i]+delta[i]/data.moveinc for i in range(len(delta))])
    data.systemCenter=data.systemCenterCopy

def testVec(vec,data): #test if a vec in cubic coords is in window
    horiztheta=math.atan(vec[1]/(vec[0]-data.windowDep)/-1)/math.pi*180
    verttheta=math.atan(vec[2]/(vec[0]-data.windowDep)/-1)/math.pi*180
    return not (abs(horiztheta)>data.horizViewTheta or( 
                abs(verttheta)>data.vertViewTheta)  or(
                vec[0]>-data.windowDep or -vec[0]>data.maxrenderdist))

def rgbStringTup(tup):
    return "#%02x%02x%02x" % tup
    
def dist(vec1,vec2):
    if len(vec1)!=len(vec2): raise Exception(
    "dist takes two n-dimensional vec's")
    return sum([(vec1[i]-vec2[i])**2 for i in range(len(vec1))])**.5
    
def init(data): 
    data.timestep=3.2E5#in seconds per frame
    data.rotateBuffer=[0,0,0]
    data.moveBuffer=[0,0,0]
    data.vecList=np.array([],dtype=float32)
    data.systemCenter=np.array([-10,0,0])
    data.systemCenterCopy=data.systemCenter
    data.numVecs=data.vecList.shape[0]
    data.filepath="picfolder" 
    data.currdate=datetime.datetime(2016,2,15,0,0,0) 
    data.defaultcurrdate=datetime.datetime(2016,2,15,0,0,0) 
    #initial conditions from arbitrary start point thanks to Solar System Viewer 
    #at https://www.fourmilab.ch/cgi-bin/Solar
    data.timestepOBJ=datetime.timedelta(seconds=data.timestep)
    data.entryString=""
    initbools(data)
    initImmutablePerams(data)
def initImmutablePerams(data):
    data.maxrenderdist=100
    data.vertViewTheta=40 #no larger than 90
    data.horizViewTheta=40
    data.psiRotInc=1
    data.thetaRotInc=1
    data.gammaRotInc=1
    data.CurrentInputs=0
    data.gravConst=6.6*10**-11
    data.windowDep=5
    data.moveinc=.1 #controls the stepsize of movement controls
def initbools(data):
    data.inEnter=False
    data.isPaused=True
    data.inIntro=True
    data.inNBod=False
    data.inSolarSim=False
    data.inHelp=False
    data.inNBodEntry=False
    data.lineTrace=True
    data.labels=False
    data.intMode=[('Euler',False),('RK4',True),('Verlet',False)] 
    data.defaultNBodVals=[['Mass [Kg]','1E20',False],
    ['Radius [Au]','.0001',False],['XPos [Au]','-10',False],
    ['YPos [Au]','1',False],['ZPos [Au]','1',False],
    ['YVel [Au/s]','2E-7',False],['XVel [Au/s]','0',False],
    ['ZVel [Au/s]','0',False],['Color','white',False]]
        
    data.firstEntry=True
def initsky():
    data.piclist={}
    skyVeclist=[]
    for file in os.listdir(data.filepath):
        if file.lower().endswith(".gif"): 
        #can only load gifs in specified directory
            picref=PhotoImage(file=data.filepath+"/"+file,width=data.width,
                              height=data.height)
            scalex=math.ceil(data.width/picref.width()/data.horizViewTheta*180)
            scaley=math.ceil(data.height/picref.height()/data.vertViewTheta*180)
            if (scalex and scaley)>1: 
                data.piclist[file]=picref.zoom(scalex,scaley) 
            else: picref.subsample(math.floor(1/scalex),math.floor(1/scaley))
    data.step=data.piclist["xneg.gif"].width()
    
def keyPressed(event, data):
    if event.keysym=='p':data.isPaused=not data.isPaused
    elif event.keysym=='h':data.inHelp=not data.inHelp
    elif event.keysym=='r' or event.keysym=='Escape' or (
        event.keysym=='c' and data.inNBod):
        resetMode()
        if event.keysym=='r':PlaceStuff()
        elif event.keysym=='Escape': 
            initbools(data)
    if data.inNBodEntry and data.isPaused:NBodInput(event,data)
    if data.inEnter:SolarEntry(event,data)
    if data.isPaused: return 
    elif event.keysym in 'qweasdzx' or event.keysym=='Up' or (
        event.keysym=='Down' or event.keysym=='Left' or event.keysym=='Right'): 
        Buffer(event,data)
    elif event.keysym=='t':data.timestep+=1E4
    elif event.keysym=='y':data.timestep-=1E4
    elif event.keysym=='i':printDebug()
    elif event.keysym=='k': data.lineTrace= not data.lineTrace
    elif event.keysym=='l': data.labels= not data.labels
def Buffer(event,data):
    if event.keysym=='q':data.rotateBuffer[0]+=1
    elif event.keysym=='e':data.rotateBuffer[0]-=1
    elif event.keysym=='w':data.rotateBuffer[1]-=1
    elif event.keysym=='a':data.rotateBuffer[2]-=1
    elif event.keysym=='s':data.rotateBuffer[1]+=1
    elif event.keysym=='d':data.rotateBuffer[2]+=1
    elif event.keysym=='z':data.moveBuffer[0]-=1
    elif event.keysym=='x':data.moveBuffer[0]+=1
    elif event.keysym=='Up':data.moveBuffer[2]+=1
    elif event.keysym=='Down':data.moveBuffer[2]-=1
    elif event.keysym=='Left':data.moveBuffer[1]-=1
    elif event.keysym=='Right':data.moveBuffer[1]+=1
def SolarEntry(event,data):
    if event.keysym=='Return': 
        data.inEnter=False
        data.defaultcurrdate=getDate(data.entryString)
        data.isPaused=False
        resetMode()
        PlaceStuff()
    elif event.keysym=='BackSpace':data.entryString=data.entryString[:-1]
    elif event.keysym=='space':data.entryString+=' '
    else: data.entryString+=event.char
    return 
def NBodInput(event,data):
        selectedList=[data.defaultNBodVals[i][2] for i in range(data.nFields)]
        if selectedList.count(True)!=1: raise Exception(
            'Must have 1 field selected!')
        selected=selectedList.index(True)
        if event.keysym=='BackSpace':
            data.defaultNBodVals[selected][1]= (
            data.defaultNBodVals[selected][1][:-1])
        elif data.firstEntry:data.defaultNBodVals[selected][1]=event.char
        else: data.defaultNBodVals[selected][1]+=event.char
        data.firstEntry=False
def resetMode():
    for planet in Planet.planets: Planet.updateFields(planet)
    for poly in Polygon.polyreg: 
        if not isinstance(poly,Planet):Polygon.updateFields(poly)

def getDate(s): #thanks python documentation!
    s+=' '
    try:
        (date,time)=tuple(s.split(' '))
        yr,mo,dy=tuple(date.split('-'))
        hr,mn,sc=tuple(time.split(':'))
        return (datetime.datetime(int(yr),int(mo),int(dy),
                int(hr),int(mn),int(sc)))
    except: 
        try: return datetime.datetime(int(yr),int(mo),int(dy))
        except: raise Exception('Invalid String')
    
def printDebug():
    print(Polygon.shapes,Polygon.shapePlaces)
    for pointIndex in Polygon.shapePlaces:
        print(data.vecList[pointIndex])
    print(data.vecList[0:10])
def mousePressed(event,data):
    if data.inIntro:
        testForStart(event,data)
    elif not data.inHelp:
        if event.x>data.intStart_x and event.y>data.intStart_y and (
            event.x<len(data.intMode)*data.rectWidth+data.intStart_x): 
            #bottom boundary is @ screen's edge
            index=(event.x-data.intStart_x)//data.rectWidth
            data.intMode[index]=(data.intMode[index][0],
                not data.intMode[index][1])
        if data.inNBod:
            testForNBodEnter(event,data)
        elif data.inSolarSim:
            if event.x>data.time_start_x and event.x<data.time_start_x+(
                data.rectWidth*data.time_scaleX) and( 
                event.y>data.time_start_y) and event.y<data.time_start_y+(
                data.rectHeight*data.time_scaleY):
                data.inEnter=True
                data.isPaused=True
def testForNBodEnter(event,data):
    if event.x>data.NBod_start_x and event.y<data.nFields*2*data.rectHeight and(
        event.x<data.NBod_start_x+data.rectWidth and event.y>data.NBod_start_y):
        index=round((event.y-data.NBod_start_y)//(2*data.rectHeight))
        data.defaultNBodVals[index][2]=not data.defaultNBodVals[index][2]
        data.inNBodEntry=True
    else:
        fields=[data.defaultNBodVals[i][0] for i in range(data.nFields)]
        XPos=data.defaultNBodVals[fields.index('XPos [Au]')][1]
        YPos=data.defaultNBodVals[fields.index('YPos [Au]')][1]
        ZPos=data.defaultNBodVals[fields.index('ZPos [Au]')][1]
        mass=data.defaultNBodVals[fields.index('Mass [Kg]')][1]
        radius=data.defaultNBodVals[fields.index('Radius [Au]')][1]
        if not data.inNBodEntry:data.p=Planet("Circle",[[float(XPos),
            float(YPos),float(ZPos)]],float(mass),float(radius),[0,0,0],
            "?","white")
        data.inNBodEntry=False
def testForStart(event,data):
        if event.x<data.rectWidth+data.solar_x and event.x>data.solar_x and(
            event.y>data.button_y and event.y<data.button_y+data.rectHeight):
            data.inIntro=False
            data.inSolarSim=True
            PlaceStuff()
        elif event.x<data.rectWidth+data.NBod_x and event.x>data.NBod_x and(
            event.y>data.button_y and event.y<data.button_y+data.rectHeight):
            data.inIntro=False
            data.inNBod=True
            PlaceStuff()    
def timerFired(data,d):#depth lets us do specific calcs on interval-calls
    if data.isPaused: return 
    oldDist,newDist=[],[]
    for planet in Planet.planets:
        if not planet.exists: continue
        oldDist+=[dist(planet.centerIsh,[0,0,0])]
    data.CurrentInputs=0
    if data.rotateBuffer!=[0,0,0]:doRotate(data.rotateBuffer[0],
        data.rotateBuffer[1],data.rotateBuffer[2],data)
    if data.moveBuffer!=[0,0,0]:doMove(data.moveBuffer[0],
        data.moveBuffer[1],data.moveBuffer[2],data)
    data.moveBuffer,data.rotateBuffer=[0,0,0],[0,0,0]
    if d==3: 
        Gravitate(data)
        if data.inSolarSim:data.currdate+=data.timestepOBJ
    for planet in Planet.planets:
        if not planet.exists: continue
        newDist+=[dist(planet.centerIsh,[0,0,0])]
    adjRad(oldDist,newDist)
    
def adjRad(list1,list2):
    passes=0
    for i in range(len(Planet.planets)):
        if not Planet.planets[i].exists: continue
        Planet.planets[i].radius*=abs(list1[passes]/list2[passes])
        passes+=1
def redrawAll(canvas, data):
    if data.inHelp: drawHelp(canvas,data)
    elif data.inIntro: drawIntro(canvas,data)
    else: 
        render(canvas, data)
        drawIntModes(canvas,data)
        if data.inNBod: drawNBodHelper(canvas,data)
        elif data.inSolarSim:
            drawTime(canvas,data) 
            if data.inEnter: drawEntry(canvas,data)

def drawIntro(canvas,data):
    introMessage="Welcome to SolarSim!\n\nPress 'h' for help at anytime!"+(
    "\n\nThis Simulation comes equipped with three modes of integration,\ni")+(
    ".e. models for calculating the effect of gravity. To simplify things,\n")+(
    "Euler integration is fast but inaccurate, RK4 is slow and accurate, but")+(
    "\ncan create physics glitches (this is because it is a nonconservative")+(
    "\ntechnique), and Verlet is fairly fast and very accurate. Better expla")+(
    "nations\nof these methods can be readily found online.\n\nPress 'Escape")+(
    "' to quit to this screen from anywhere.\nHave Fun!")
    intro_x,intro_y,Ifont=data.width//2,data.height//8,"Helvetica 14"
    solar_x=data.width//4
    button_y=3*data.height//4
    rectCol='grey'
    rectWidth,rectHeight=80,20
    NBod_x=data.width//4*3
    canvas.create_text(intro_x,intro_y,anchor=N,text=introMessage,font=Ifont)    
    canvas.create_rectangle(solar_x,button_y,rectWidth+solar_x,
        rectHeight+button_y,fill=rectCol)
    canvas.create_text(solar_x+rectWidth/2,rectHeight/2+button_y,
        text="SolarSim")
    canvas.create_rectangle(NBod_x,button_y,rectWidth+NBod_x,
        rectHeight+button_y,fill=rectCol)
    canvas.create_text(NBod_x+rectWidth/2,rectHeight/2+button_y,text="N-Body")
    data.rectWidth,data.rectHeight,data.rectCol=rectWidth,rectHeight,rectCol
    data.solar_x,data.button_y=solar_x,button_y
    data.NBod_x=NBod_x
    
def drawHelp(canvas,data):
    help_x=data.width//2
    help_y=data.height/8
    control_x=data.width//2
    control_y=4*data.height//9
    helpfont=('Arial',14)
    helpmsg="SolarSim has two operational modes: An N-Body gravitation tester"+(
    ", \nand a Solar simulator. Controls are (mostly) shared between each, a")+( 
    "nd are as follows:")
    controlmsg=" 'a': yaw left \n 'd': yaw right \n 'e': roll C.W. \n" +(
        "'i': debug info. \n 'p': pause/unpause \n 'q': roll C.C.W \n 'r': ")+( 
        "restart game \n's': pitch down \n 't': inc. timestep (by 1000 [S])")+(
        "\n 'w': pitch up \n'x': zoom in \n 'y': dec. timestep (by 1000 [S])")+(
        "\n 'z': zoom out \n 'Up': move up \n 'Down': move down \n'Left': ")+( 
        "move left \n 'Right': move right \n 'c': clear planets (NBod only)")+(
        " \n 'k': enable/disable trajectory tracing \n 'l': show/hide labels")+(
        "(SolarSim only)")
    canvas.create_text(help_x,help_y,text=helpmsg,font=helpfont)
    canvas.create_text(control_x,control_y,text=controlmsg,font=helpfont)
    userguide=(
    "Some notes for use: \n-In the SolarSim, dates are given in Year-Month-Da"+(
    "y Hour:Minute:Second format\n-To create a new planet in the N-body Sim.")+( 
    ",select the field you want to\nedit and double-click anywhere to place")+( 
    " when finished editing\n-To create a custom, non-random scenario,\njust")+( 
    "press 'c' in the N-Body Sim. and design new planets"))
    canvas.create_text(control_x,control_y*2,text=userguide,font=helpfont)
def drawIntModes(canvas,data):
    rWidth=data.rectWidth
    rHeight=data.rectHeight
    modeCopy=data.intMode
    num=len(modeCopy)
    start_x=0
    start_y=data.height-rHeight
    fillcols=['grey' if modeCopy[i][1] else 'white' for i in range(num)]
    for i in range(len(modeCopy)):
        canvas.create_rectangle(start_x+rWidth*i,start_y,start_x+rWidth*(i+1),
            start_y+rHeight, fill=fillcols[i])
        canvas.create_text(start_x+rWidth*(i+.5),start_y+rHeight/2,
            text=modeCopy[i][0])
    data.intStart_x=start_x
    data.intStart_y=start_y
def drawNBodHelper(canvas,data):
    localvals=copy.deepcopy(data.defaultNBodVals)
    start_x=data.width-data.rectWidth
    start_y=0
    nFields=len(localvals)
    rectcol=['grey' if localvals[i][2] else 'white' for i in range(nFields)]
    rectHeight=data.rectHeight*2
    for i in range(nFields):
        canvas.create_rectangle(start_x,start_y+rectHeight*i,
            start_x+data.rectWidth,start_y+rectHeight*(i+1),fill=rectcol[i])
        canvas.create_text(start_x+data.rectWidth*.5,start_y+rectHeight*i,
            anchor=N,text=localvals[i][0])
        canvas.create_text(start_x+data.rectWidth*.5,start_y+rectHeight*(i+1),
            anchor=S,text=localvals[i][1])
    data.NBod_start_x,data.NBod_start_y,data.nFields=start_x,start_y,nFields

def drawTime(canvas,data):
    currdate=str(data.currdate)
    start_x,scaleX=0,3
    start_y,scaleY=0,2
    rectcol='white'
    rectHeight=data.rectHeight*scaleY
    rectWidth=data.rectWidth*scaleX
    step=data.timestep
    stepText='Stepping at %.2e Sec/Tick' % step
    anchorYPos=start_y+rectHeight/2
    anchorXPos=start_x+rectWidth/2
    canvas.create_rectangle(start_x,start_y,start_x+rectWidth,
        start_y+rectHeight,fill=rectcol)
    canvas.create_text(anchorXPos,anchorYPos,anchor=N,text=stepText)
    canvas.create_text(anchorXPos,anchorYPos,anchor=S,text=currdate)
    data.time_start_x,data.time_start_y=start_x,start_y
    data.time_scaleX,data.time_scaleY=scaleX,scaleY
    
def drawEntry(canvas,data):
    scale=3
    entryhelp="Enter a valid date of the form" +(
    "<Year-Month-Day Hour:Minute:Second>\n")+(
    "                            (A date is required, time is optional)")
    canvas.create_rectangle(data.width//2-scale*data.rectWidth,
        data.height//2-2*data.rectHeight,data.width//2+scale*data.rectWidth,
        data.height//2+2*data.rectHeight,fill='white')
    yoffset=10
    canvas.create_text(data.width//2,-yoffset+data.height//2,anchor=S,
        text=entryhelp)
    canvas.create_text(data.width//2,data.height//2+yoffset,anchor=N,
        text=data.entryString)

def render(canvas,data):
    canvas.create_image(0,0,image=data.piclist["xneg.gif"])
    linecol='white'
    passes=0
    for poly in Polygon.shapePlaces.keys():
        if not Polygon.polyreg[poly].exists:continue
        reflist=Polygon.shapePlaces[poly]
        if Polygon.shapes[poly]=="Circle":
            center=data.vecList[reflist[0]]
            if not testVec(center,data): continue
            drawCircle(passes,center,poly,data,canvas)
            passes+=1
        else:
            for i in range(-1,len(reflist)-1):
                vec=data.vecList[reflist[i]]
                nextvec=data.vecList[reflist[i+1]]
                if not testVec(vec,data):
                    continue
                scaledvec=scaleVec(vec,data)
                adjx=data.width/2 -scaledvec[1]#adjusted x and y from y and 
                adjy=data.height/2+scaledvec[2]#z coords in cubic
                nextx=data.width/2-scaleVec(nextvec,data)[1]
                nexty=data.height/2+scaleVec(nextvec,data)[2]
                canvas.create_line(adjx,adjy,nextx,nexty,fill=linecol)
    cx,cy,r=data.width//2,data.height//2,3 #this is the cursor
    canvas.create_oval(cx-r,cy-r,cx+r,cy+r,fill='white',width=1,
                       outline='black')
def drawCircle(passes,center,poly,data,canvas): #Planets NOT drawn to scale
    keystart_x=data.width*4/5
    keystart_y=data.height*9/10
    keyscale=40
    expscale=.2
    multscale=40
    rings=8
    adjx=data.width/2-scaleVec(center,data)[1]
    adjy=data.height/2+scaleVec(center,data)[2]
    r=multscale*Polygon.polyreg[poly].radius**expscale
    color=Polygon.polyreg[poly].color
    if isinstance(Polygon.polyreg[poly],Planet) and( 
        Polygon.polyreg[poly].name=='Saturn'):
        canvas.create_oval(adjx-r-rings,adjy-r- rings,adjx+r+rings,
            adjy+r+rings,fill=color)
        canvas.create_oval(adjx-r-rings+2,adjy-r+2-rings,adjx+r+rings-2,
            adjy+r+rings-2,fill='black')
    canvas.create_oval(adjx-r,adjy-r,adjx+r,adjy+r,fill=color)
    if isinstance(Polygon.polyreg[poly],Planet) and data.inSolarSim:
        if data.labels:canvas.create_text(adjx+r,adjy+r,anchor=NW,
            text=Polygon.polyreg[poly].name,fill='white')
        canvas.create_oval(keystart_x-r,keystart_y-passes*keyscale-r,
            keystart_x+r,keystart_y-passes*keyscale+r,fill=color)
        canvas.create_text(keystart_x+keyscale,keystart_y-passes*keyscale,
            anchor=W,text=Polygon.polyreg[poly].name, fill='white')
    
def runrender(width=800,height=800):#This is thanks to Dr. Kosbie's Tkinter Tech
    def redrawAllWrapper(canvas, data):
        canvas.delete(ALL)
        redrawAll(canvas, data)
        canvas.update()    

    def mousePressedWrapper(event, canvas, data):
        mousePressed(event, data)
        data.CurrentInputs+=1

    def keyPressedWrapper(event, canvas, data):
        keyPressed(event, data)
        data.CurrentInputs+=1  

    def timerFiredWrapper(canvas, data,d=0):
        d+=1
        d%=6
        timerFired(data,d)
        redrawAllWrapper(canvas, data)
        canvas.after(data.timerDelay, timerFiredWrapper, canvas, data,d)
    Planet.clearFields()
    root = Toplevel()
    data.width = 800
    data.height = 800
    data.timerDelay = 120 # milliseconds 
    initsky()
    canvas = Canvas(root, width=data.width, height=data.height)
    canvas.pack()
    root.bind("<Button-1>", lambda event:
                            mousePressedWrapper(event, canvas, data))
    root.bind("<Key>", lambda event:
                            keyPressedWrapper(event, canvas, data))
    timerFiredWrapper(canvas, data)
    root.mainloop() 
class Struct(object): 
    pass    
data = Struct()
init(data)
 
class Polygon(object):
    Polycount=0
    shapes={}
    shapePlaces={}
    polyreg=[]
    def __init__(self,shape,pointmat,color,radius=0):
        self.color,self.radius,self.exists,self.shape=color,radius,True,shape
        self.ID=Polygon.Polycount
        initHelper(pointmat)
        self.centerIsh=np.array([0,0,0],dtype=float64)
        vecListStarter=len(data.vecList)
        nGon=len(pointmat)
        Polygon.shapes[self.ID]=shape
        Polygon.polyreg.append(self)
        for pointIndex in range(nGon):
            if pointIndex!=0:
                Polygon.shapePlaces[self.ID]+=[vecListStarter+pointIndex]
            else: Polygon.shapePlaces[self.ID]=[vecListStarter+pointIndex]
            self.centerIsh[0]+=pointmat[pointIndex][0]/nGon
            self.centerIsh[1]+=pointmat[pointIndex][1]/nGon
            self.centerIsh[2]+=pointmat[pointIndex][2]/nGon
        if data.vecList!=[]:
            data.vecList=np.vstack((data.vecList,pointmat)) 
        else: data.vecList=np.array(pointmat,dtype=float32)
        data.numVecs=data.vecList.shape[0]

    def updateFields(self): #cheating!
        self.exists=False
    def __repr__(self):
        return "%s at %s" % (self.shape, self.points)
def initHelper(pointmat):
    try: len(pointmat)
    except: raise Exception("Points must given in a sliceable format")
    Polygon.Polycount+=1
    
def PlaceStuff(): #style is pretty bad here, but I think this function is actually
    #as simple as it could be--could split 80+char lines, but they are all
    #info for the same planet so this would detract from clarity in my opinion
    SunMass,SunRad=2E30,.047
    EarthMass,EarthRad,EarthStartTheta=6E24,4.3E-5,100
    toRad=180/math.pi
    toDeg=360
    if data.inSolarSim: #planetary info from 
        #http://www.windows2universe.org/our_solar_system/planets_table.html
        #in (estimated) standard polar-form-degrees from 2016-01-01 0:00:00
        SecPerYr=3.154e+7
        timechange=(datetime.datetime(2016,1,1)-data.defaultcurrdate).total_seconds()
        timechange*=SecPerYr
        MercuryYr,MercuryStartTheta,MercuryAU,MercuryMass,MercuryRad=.24,30,.39,.055*6E24,1.63E-5
        VenusYr,VenusStartTheta,VenusAU,VenusMass,VenusRad=.62,185,.72,.815*6E24,4.04E-5
        MarsYr,MarsStartTheta,MarsAU,MarsMass,MarsRad=1.88,170,1.52,.107*6E24,2.27E-5
        JupiterYr,JupiterStartTheta,JupiterAU,JupiterMass,JupiterRad=11.86,160,5.2,318*6E24,0.00048
        SaturnYr,SaturnStartTheta,SaturnAU,SaturnMass,SaturnRad=29.46,240,9.54,95*6E24,0.0004
        UranusYr,UranusStartTheta,UranusAU,UranusMass,UranusRad=84.01,20,19.18,15*6E24,0.000171
        NeptuneYr,NeptuneStartTheta,NeptuneAU,NeptuneMass,NeptuneRad=164.8,330,30.06,17*6E24,0.00016
        earthcol=rgbStringTup((0,63,255))
        MercuryTheta=(MercuryStartTheta+timechange/MercuryYr*toDeg)/toRad
        EarthTheta=(EarthStartTheta+timechange*toDeg)/toRad
        VenusTheta=(VenusStartTheta+timechange/VenusYr*toDeg)/toRad
        MarsTheta=(MarsStartTheta+timechange/MarsYr*toDeg)/toRad
        JupiterTheta=(JupiterStartTheta+timechange/JupiterYr*toDeg)/toRad
        SaturnTheta=(SaturnStartTheta+timechange/SaturnYr*toDeg)/toRad
        UranusTheta=(UranusStartTheta+timechange/UranusYr*toDeg)/toRad
        NeptuneTheta=(NeptuneStartTheta+timechange/NeptuneYr*toDeg)/toRad
        EarthX,EarthY=math.cos(EarthTheta),math.sin(EarthTheta)
        MarsX,MarsY=MarsAU*math.cos(MarsTheta),MarsAU*math.sin(MarsTheta)
        JupiterX,JupiterY=JupiterAU*math.cos(JupiterTheta),JupiterAU*math.sin(JupiterTheta)
        SaturnX,SaturnY=SaturnAU*math.cos(SaturnTheta),SaturnAU*math.sin(SaturnTheta)
        UranusX,UranusY=UranusAU*math.cos(UranusTheta),UranusAU*math.sin(UranusTheta)
        NeptuneX,NeptuneY=NeptuneAU*math.cos(NeptuneTheta),NeptuneAU*math.sin(NeptuneTheta)
        MercuryX,MercuryY=MercuryAU*math.cos(MercuryTheta),MercuryAU*math.sin(MercuryTheta)
        VenusX,VenusY=VenusAU*math.cos(VenusTheta),VenusAU*math.sin(VenusTheta)
        p1=Planet("Circle",[[-10,EarthX,EarthY]],EarthMass,EarthRad,
            getVel([0,EarthX,EarthY],1),"Earth",earthcol)
        p2=Planet("Circle",[data.systemCenterCopy],SunMass,
            SunRad,[0,0,0],"Sun","orange")
        p3=Planet("Circle",[[-10,VenusX,VenusY]],VenusMass,VenusRad,
            getVel([0,VenusX,VenusY],VenusYr),"Venus","green")
        p4=Planet("Circle",[[-10,MercuryX,MercuryY]],MercuryMass,MercuryRad,
            getVel([0,MercuryX,MercuryY],MercuryYr),"Mercury","yellow")
        p5=Planet("Circle",[[-10,MarsX,MarsY]],MarsMass,MarsRad,
            getVel([0,MarsX,MarsY],MarsYr),"Mars","red")
        p6=Planet("Circle",[[-10,JupiterX,JupiterY]],JupiterMass,JupiterRad,
            getVel([0,JupiterX,JupiterY],JupiterYr),"Jupiter","orange")
        p7=Planet("Circle",[[-10,SaturnX,SaturnY]],SaturnMass,SaturnRad,
            getVel([0,SaturnX,SaturnY],SaturnYr),"Saturn","yellow")
        p8=Planet("Circle",[[-10,NeptuneX,NeptuneY]],NeptuneMass,NeptuneRad,
            getVel([0,NeptuneX,NeptuneY],NeptuneYr),"Neptune","dark blue")
        p9=Planet("Circle",[[-10,UranusX,UranusY]],UranusMass,UranusRad,
            getVel([0,UranusX,UranusY],UranusYr),"Uranus","blue")
    elif data.inNBod:
        N=random.randint(2,5)
        planetlist=[] 
        Sun=Planet("Circle",[data.systemCenterCopy],SunMass,SunRad,[0,0,0],
                    "Sun","orange")
        for currentPlanet in range(N): 
            mass=random.randint(1,100)*10**random.randint(15,30)
            rad=random.randint(0,100000)/100000000
            x,y,z=-10,random.randint(-50,50)/10,random.randint(-50,50)/10
            vx,vy,vz=0,random.randint(-10,10)/1E8,random.randint(-50,50)/1E8
            p=Planet("Circle",[[x,y,z]],mass,rad,[vx,vy,vz])
def getVel(Rad,Tau):
    SecPerYr=3.154e+7
    omega=2*math.pi/(Tau*SecPerYr)
    return [Rad[0]*omega,-Rad[2]*omega,Rad[1]*omega] # all the planets rotate c.w.
class Planet(Polygon):
    planets=[]
    def __init__(self,shape,CoM,mass,radius,vel,name='Unknown',fillcol='White'):
        super().__init__(shape,CoM,fillcol,radius)
        self.mass=mass
        self.vel=np.array(vel,dtype=float32)
        Planet.planets.append(self)
        self.forces=[]
        self.name=name
        self.started=True
        self.lastPos=np.array([0,0,0]) #resets
    def __repr__(self):
        return "%s (%s) at %s" % (self.name,self.shape,self.centerIsh)
    @staticmethod
    def clearFields():
        Polygon.shapes={}
        Polygon.shapePlaces={}
        Polygon.Polycount=0
        planets=[]
        polyreg=[]
    def checkForCollision(self,other): #using super lame bounding spheres, 
    #works fine for this project
        if not isinstance(other,Planet): 
            raise Exception("NonPlanet-Collision")
        if not (self.exists and other.exists): return False
        nSteps=4
        collided=False
        starter=1 if self.started else 0
        path1=self.centerIsh-self.lastPos
        path2=other.centerIsh-other.lastPos
        for step1 in range(starter,nSteps+1):
            pos1=self.lastPos+step1/nSteps*path1
            for step2 in range(starter,nSteps+1):
                pos2=other.lastPos+step2/nSteps*path2
                if dist(pos1,pos2)<other.radius+self.radius:
                    collided=True
                    break
        return (dist(self.centerIsh,other.centerIsh)<=(self.radius+(
                other.radius)) or collided)
    def __eq__(self,other):
        if not isinstance(other,Planet): 
            raise Exception("Can't compare planet and non-planet")
        return self.name==other.name and self.ID==other.ID
    def sumForces(self):
        localforce=np.array([0,0,0],dtype=float64)
        if len(self.forces)==0: 
            self.forces=[(0,localforce)]
            return 
        if self.name!='Sun':
            for force in self.forces:
                localforce=np.ndarray.astype(localforce+force[0]*force[1],
    dtype=float64,casting='unsafe')
        localmag=mag(localforce) if mag(localforce)!=0 else 1
        self.forces=[(localmag,localforce/localmag)]
    def updateFields(deleted): #currently just keeping veclist entries 
    #(not really deleting)
        super().updateFields()
        deleted.mass=1
        deleted.radius=.00000000001
        deleted.color='red'
        deleted.reflist=[]
    @staticmethod
    def combineBodies():
        checked=[]
        for planet1 in Planet.planets:
            if not planet1.exists: continue
            for planet2 in Planet.planets:
                if (planet2,planet1) in checked or planet1==planet2:
                    continue
                checked.append((planet1,planet2))
                if Planet.checkForCollision(planet1,planet2):
                    d1=planet1.mass/(4*math.pi*planet1.radius**3)
                    d2=planet2.mass/(4*math.pi*planet2.radius**3)
                    vFinal=(planet1.vel*planet1.mass+planet2.vel*planet2.mass)/(
                            planet1.mass+planet2.mass)
                    if planet1.mass>planet2.mass:
                        planet1.vel=vFinal
                        planet1.mass+=planet2.mass
                        planet1.radius=(planet1.mass/d1/4/math.pi)**(1/3)
                        Planet.updateFields(planet2)
                    else:
                        planet2.vel=vFinal
                        planet2.mass+=planet1.mass
                        planet2.radius=(planet2.mass/d2/4/math.pi)**(1/3)
                        Planet.updateFields(planet1)
def mag(vec):
    return sum([vec[i]**2 for i in range(len(vec))])**.5
def Force(m1,m2,r):
    return m1*m2*data.gravConst/r/r 
def updateForces(data):
    minForce=1E3 #[N]
    maxForce=1E30
    refSum=0
    MPerAU=1.496e+11
    for planet1 in Planet.planets:
        if not planet1.exists:continue
        if len(Planet.planets)<2:
            break
        for ref in Polygon.shapePlaces[planet1.ID]:
            refSum+=data.vecList[ref]
        planet1.centerIsh=refSum/len(Polygon.shapePlaces[planet1.ID])
        refSum=0
        for planet2 in Planet.planets:
            if planet1==planet2 or not planet2.exists:
                continue
            Rad=dist(planet1.centerIsh,planet2.centerIsh) #[Au]
            forcemag=Force(planet1.mass,planet2.mass,Rad*MPerAU) #[N]
            forcedir=(planet2.centerIsh-planet1.centerIsh)/Rad
            if forcemag>minForce:planet1.forces.append((forcemag,forcedir))
            elif forcemag>maxForce: planet1.forces.append((maxForce, forcedir))

def Gravitate(data):
    updateForces(data)
    MPerAU=1.496e+11
    for planet in Planet.planets:
        if not planet.exists: continue
        Planet.sumForces(planet) 
        templastPos=planet.centerIsh
        A=planet.forces[0][0]*planet.forces[0][1]/planet.mass/MPerAU #[Au/s**2]
        planet.forces=[]
        if data.intMode[0][1]: nextState=(planet.centerIsh+(
                        planet.vel*data.timestep),planet.vel+A*data.timestep)
        elif data.intMode[1][1]:nextState=RK4(planet.centerIsh,
                        planet.vel,A,data.timestep,planet)
        elif data.intMode[2][1]:nextState=Verlet(planet,A,data)
        else: raise Exception("No Integration Mode Selected!!!")
        planet.centerIsh,planet.vel=nextState
        planet.lastPos=templastPos
        Planet.combineBodies()
        if data.lineTrace: Polygon("line",
            [planet.lastPos,planet.centerIsh],'white')
        for ref in Polygon.shapePlaces[planet.ID]:
            data.vecList[ref]=planet.centerIsh

def Verlet(planet,A,data): #explanation from wikipedia under 'Verlet'
    if not planet.started: 
        nextR=(2*planet.centerIsh-planet.lastPos+A*(data.timestep**2))
    else: 
        nextR=planet.centerIsh+planet.vel*data.timestep+A*(data.timestep**2)/2
        planet.started=False
    nextV=(nextR-planet.centerIsh)/data.timestep
    nextState=(nextR,nextV)
    return nextState
def RK4(r,v,a,step,planet): #aided by this paper's explanation: 
    #http://spiff.rit.edu/richmond/nbody/OrbitRungeKutta4.pdf
    def getAccelfromNewPos(planet1,pos):
        MPerAU=1.496e+11
        localforce=np.array([0,0,0],dtype=float32)
        for planet2 in Planet.planets:
            if planet1==planet2 or not planet2.exists:
                continue
            Rad=dist(planet1.centerIsh+pos,planet2.centerIsh)
            localforce=np.ndarray.astype(localforce+(
                Force(planet1.mass,planet2.mass,Rad*MPerAU)*(
                planet2.centerIsh-planet1.centerIsh)/Rad),
    dtype=float64,casting='unsafe')
        return localforce/planet1.mass/MPerAU
    k1r=v
    k1v=a
    k2r=v*step/2*k1v
    k2v=getAccelfromNewPos(planet,k1r*step/2)
    k3r=v*step/2*k2v
    k3v=getAccelfromNewPos(planet,k2r*step/2)
    k4r=v*step*k3v
    k4v=getAccelfromNewPos(planet,k3r*step)
    nextV=v+step/6*(k1v+2*k2v+2*k3v+k4v)
    nextR=r+step/6*(k1r+2*k2r+2*k3r+k4r) 
    return (nextR,nextV)
runrender()






