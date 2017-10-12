from sympy.geometry import *
from sympy import pi,symbols,atan2
from sympy.plotting import plot_implicit
from matplotlib import pyplot as plt
import random

def plot_VD_CH(VD,CH,points):
	# print(CH)
	plt.clf()
	for p in points:
		plt.plot([p.x],[p.y],marker='o',markersize=3,color='blue')
		plt.annotate(	"({0},{1})".format(p.x,p.y),
						xy=(p.x, p.y),
						xytext=(5, 5),
						textcoords='offset points',
						# ha='right',
						# va='bottom',
						)

	plt.plot([p.x for p in CH]+[CH[0].x,], [p.y for p in CH]+[CH[0].y,], 'r--')
	if(len(points)>1):
		Xmin=min(*[p.x for p in CH])
		Ymin=min(*[p.y for p in CH])
		Xmax=max(*[p.x for p in CH])
		Ymax=max(*[p.y for p in CH])
	else:
		Xmin=points[0].x
		Ymin=points[0].y
		Xmax=points[0].x
		Ymax=points[0].y
	for v in VD.values():
		if not isinstance(v,Point2D):
			Xmin = min(Xmin,v.p1.x,v.p2.x)
			Ymin = min(Ymin,v.p1.y,v.p2.y)
			Xmax = max(Xmax,v.p1.x,v.p2.x)
			Ymax = max(Ymax,v.p1.y,v.p2.y)
	# print (Xmin,Ymin,Xmax,Ymax)
	for k,v in VD.items():
		# print (k,"\t",v)
		if isinstance(v,Segment2D):
			plt.plot([v.p1.x,v.p2.x],[v.p1.y,v.p2.y],color = 'black')
		elif isinstance(v,Ray2D):
			a,b,c = Line2D(v).coefficients
			if b:
				if v.direction.x<0:
					y=(-c-a*Xmin)/b
					if y<=Ymax:
						if y>=Ymin:
							plt.plot([v.source.x,Xmin],[v.source.y,(-c-a*Xmin)/b],color = 'black')
						else:
							plt.plot([v.source.x,(-c-b*Ymin)/a],[v.source.y,Ymin],color = 'black')
					else:
						plt.plot([v.source.x,(-c-b*Ymax)/a],[v.source.y,Ymax],color = 'black')
				else:
					y=(-c-a*Xmax)/b
					if y<=Ymax:
						if y>=Ymin:
							plt.plot([v.source.x,Xmax],[v.source.y,(-c-a*Xmax)/b],color = 'black')
						else:
							plt.plot([v.source.x,(-c-b*Ymin)/a],[v.source.y,Ymin],color = 'black')
					else:
						plt.plot([v.source.x,(-c-b*Ymax)/a],[v.source.y,Ymax],color = 'black')
			else:
				if v.direction.y<0:
					plt.plot([v.source.x,v.source.x],[v.source.y,Ymin],color = 'black')
				else:
					plt.plot([v.source.x,v.source.x],[v.source.y,Ymax],color = 'black')
		elif isinstance(v,Line2D):
			a,b,c = v.coefficients
			if b:
				x1,y1 = Xmin,(-c-a*Xmin)/b
				if y1>Ymax:
					x1,y1 = (-c-b*Ymax)/a,Ymax
				elif y1<Ymin:
					x1,y1 = (-c-b*Ymin)/a,Ymin
				x2,y2 = Xmax,(-c-a*Xmax)/b
				if y2>Ymax:
					x2,y2 = (-c-b*Ymax)/a,Ymax
				elif y2<Ymin:
					x2,y2 = (-c-b*Ymin)/a,Ymin
				plt.plot([x1,x2],[y1,y2],color = 'black')
			else:
				plt.plot([v.p1.x,v.p1.x],[Ymax,Ymin],color = 'black')
		elif isinstance(v,Point2D):
			plt.plot([v.x,],[v.y,],'o')

	plt.xlim(-20,20)
	plt.ylim(-20,20)
	plt.show()


def indexOfRightMostPoint(points):
	a=[each.x for each in points]
	return a.index(max(a))

def indexOfLeftMostPoint(points):
	a=[each.x for each in points]
	return a.index(min(a))

def getUpperTangent(leftCH,rightCH):
	llen = len(leftCH)
	rlen = len(rightCH)
	il = indexOfRightMostPoint(leftCH)
	ir = indexOfLeftMostPoint(rightCH)

	flag = True
	while flag:
		flag = False

		temp = leftCH
		while temp:
			tangent = Segment(leftCH[il],rightCH[ir])
			a,b,c = Line(tangent).coefficients
			temp = [ p for p in temp if ((b*(a*p.x + b*p.y +c)) > 0 if b!=0 
																	else ((p.y > leftCH[il].y) if (p.x == leftCH[il].x)
																								else True
																	)
										)
					]
			if temp:
				il = (il-1) % llen
				flag = True

		temp = rightCH
		while temp:
			tangent = Segment(leftCH[il],rightCH[ir])
			a,b,c = Line(tangent).coefficients
			temp = [p for p in temp if ((b*(a*p.x + b*p.y +c)) > 0 if b!=0 
																	else ( (p.y < rightCH[ir].y) if (p.x == rightCH[ir].x)
																									else False
																	)
									)
					]
			if temp:
				ir = (ir+1) %rlen
				flag = True
			
	return tangent

def getlowerTangent(leftCH,rightCH):
	llen = len(leftCH)
	rlen = len(rightCH)
	il = indexOfRightMostPoint(leftCH)
	ir = indexOfLeftMostPoint(rightCH)

	flag = True
	while flag:
		flag = False

		temp = leftCH
		while temp:
			tangent = Segment(leftCH[il],rightCH[ir])
			a,b,c = Line(tangent).coefficients
			temp = [ p for p in temp if ((b*(a*p.x + b*p.y +c)) < 0 if b!=0 
																	else ((p.y > leftCH[il].y) if (p.x == leftCH[il].x)
																								else False
																	)
										)
					]
			if temp:
				il = (il+1) % llen
				flag = True

		temp = rightCH
		while temp:
			tangent = Segment(leftCH[il],rightCH[ir])
			a,b,c = Line(tangent).coefficients
			temp = [p for p in temp if ((b*(a*p.x + b*p.y +c)) < 0 if b!=0 
																	else ( (p.y < rightCH[ir].y) if (p.x == rightCH[ir].x)
																									else True
																	)
									)
					]
			if temp:
				ir = (ir-1) %rlen
				flag = True
			
	return tangent

def VDKey_intersection_at_highest_point(line,VD):
	result=None
	pointofIntersection=None
	for key,value in VD.items():
		intersection = value.intersection(line)
		# print ("hp\t",value ,"\t", line,"\t",intersection)
		if intersection  :
			if isinstance(intersection[0], Point2D):
				topMostPoint=intersection[0]
			elif isinstance(intersection[0] , Ray2D): 
				topMostPoint=intersection[0].source
			elif isinstance(intersection[0] , Segment2D):
				topMostPoint=intersection[0].p1 if (intersection[0].p1.y>intersection[0].p2.y or (intersection[0].p1.y==intersection[0].p2.y and intersection[0].p1.x<intersection[0].p2.x)) else intersection[0].p2
			else:
				continue
			if (result == None or (topMostPoint.y>pointofIntersection.y or (topMostPoint.x<pointofIntersection.x and (topMostPoint.y==pointofIntersection.y)))):
				result =key
				pointofIntersection = topMostPoint
	return result,pointofIntersection

def lie_left(BisectingLines,point):
	if point.y>= BisectingLines[0][1].source.y:
		start=0
	elif point.y <=BisectingLines[-1][1].source.y:
		start=-1
	else:
		start=1
		end=len(BisectingLines)-1
		while(start<end):
			mid=int((start+end)/2)
			Ymax=max(BisectingLines[mid][1].p1.y,BisectingLines[mid][1].p2.y)
			Ymin=min(BisectingLines[mid][1].p1.y,BisectingLines[mid][1].p2.y)
			if(point.y<=Ymax and point.y>=Ymin):
				start=mid
				break
			elif (point.y > Ymax):
				end=mid
			else:
				start=mid+1
	a,b,c = Line2D(BisectingLines[start][1]).coefficients
	return point.x <= (min(BisectingLines[start][1].p1.x,BisectingLines[start][1].p2.x) if a==0 else ((-c-b*point.y)/a))

def mergeVDUtil(leftVD,rightVD,BisectingLines):
	# print ("left:",leftVD)
	# print ("right:",rightVD)
	VD=dict((each[0],each[1]) for each in BisectingLines)
	flag = (len(BisectingLines)==1)
	for eachKey,eachValue in leftVD.items():
		if flag or (lie_left(BisectingLines,eachValue.p1) and lie_left(BisectingLines,eachValue.p2)):
			VD[eachKey]=eachValue
	for eachKey,eachValue in rightVD.items():
		if flag or (not lie_left(BisectingLines,eachValue.p1)) or (not lie_left(BisectingLines,eachValue.p2)):
			VD[eachKey]=eachValue
	return VD

def voronoiLinesUtil(points):
	npoints = len(points)

	if npoints == 1:
		return {},[points[0],]
	if points[0].x==points[-1].x    or    npoints<=2:
		return {(points[i],points[i+1]) : Segment(points[i],points[i+1]).perpendicular_bisector() for i in range(npoints-1)},[points[0],points[-1]]

	llen=int(npoints/2)
	rlen=npoints-llen
	# print(points)
	leftVD,leftCH = voronoiLinesUtil(points[:llen])
	# plot_VD_CH(leftVD,leftCH,points[:llen])
	rightVD,rightCH = voronoiLinesUtil(points[llen:])
	# plot_VD_CH(rightVD,rightCH,points[llen:])
	# print("end")
	tangent = getUpperTangent(leftCH,rightCH)
	# print ("upper tangnet\t",tangent)
	upperTangent = tangent
	# note: I have sorted points on the bases of there x and y coordinates so for any tangent p1 lies in leftCH and p2 lies in rightCH
	
	bisectingLine = tangent.perpendicular_bisector()
	if atan2(bisectingLine.direction.y,bisectingLine.direction.x) > 0:
		bisectingLine=Line2D(bisectingLine.p2,bisectingLine.p1)
	# print (leftVD)
	# print (rightVD)
	lIVDKey,lPoI = VDKey_intersection_at_highest_point(bisectingLine,leftVD)
	rIVDKey,rPoI = VDKey_intersection_at_highest_point(bisectingLine,rightVD)
	# print ("poi")
	# print (lIVDKey,lPoI)
	# print (rIVDKey,rPoI)

	prePoI=None
	lines_from_PoI_left={}
	lines_from_PoI_right={}
	BisectingLines=[]

	while lPoI or rPoI:
		if lPoI and rPoI:
			if lPoI.y>rPoI.y or (lPoI.y==rPoI.y and lPoI.x<=rPoI.x):
				PoI,IVDKey,IVDLine=lPoI,lIVDKey,leftVD[lIVDKey]
				del leftVD[lIVDKey]
			else:
				PoI,IVDKey,IVDLine=rPoI,rIVDKey,rightVD[rIVDKey]
				del rightVD[rIVDKey]
		elif lPoI:
			PoI,IVDKey,IVDLine=lPoI,lIVDKey,leftVD[lIVDKey]
			del leftVD[lIVDKey]
		else:
			PoI,IVDKey,IVDLine=rPoI,rIVDKey,rightVD[rIVDKey]
			del rightVD[rIVDKey]

		if(PoI!=prePoI):
			prePoI=PoI
			for eachKey,eachValue in lines_from_PoI_left.items():
				leftVD[eachKey]=eachValue
			for eachKey,eachValue in lines_from_PoI_right.items():
				rightVD[eachKey]=eachValue
			lines_from_PoI_left.clear()
			lines_from_PoI_right.clear()

		# bisectingline
		if isinstance(bisectingLine,Line2D):
			BisectingLines.append(((tangent.p1,tangent.p2) , Ray(PoI,PoI-bisectingLine.direction)))
		elif isinstance(Segment(PoI,bisectingLine.source),Segment2D):
			BisectingLines.append(((tangent.p1,tangent.p2) , Segment(PoI,bisectingLine.source)))


		# creating VD
		l1=Ray(PoI,PoI+IVDLine.direction)
		l1_in_left = (atan2(l1.direction.y,l1.direction.x)-atan2(bisectingLine.direction.y,bisectingLine.direction.x))
		while l1_in_left > pi:	l1_in_left -= 2*pi
		while l1_in_left <= -pi:	l1_in_left += 2*pi
		l1_in_left = (l1_in_left <= 0)

		if (l1_in_left) == (PoI == lPoI):
			if isinstance(IVDLine,Segment2D):
				l1=Segment(PoI,IVDLine.p2)
		else:
			if isinstance(IVDLine,Line2D):
				l1=Ray(PoI,PoI-IVDLine.direction)
			elif isinstance(IVDLine,Ray2D):
				l1=Segment(PoI,IVDLine.source)
			else:
				l1=Segment(PoI,IVDLine.p1)

		# print(PoI,lPoI,rPoI)
		# print(IVDKey,tangent)
		# print("\n")
		if PoI == lPoI:
			if not isinstance(l1,Point2D):
				lines_from_PoI_left[IVDKey]=l1
			tangent = Segment(IVDKey[0] if tangent.p1 == IVDKey[1] else IVDKey[1],tangent.p2)
		else:
			if not isinstance(l1,Point2D):
				lines_from_PoI_right[IVDKey]=l1
			tangent = Segment(tangent.p1,IVDKey[0] if tangent.p2 == IVDKey[1] else IVDKey[1])
		# print(IVDKey,tangent)
		# print("\n")
		bisectingLine = tangent.perpendicular_bisector()
		if atan2(bisectingLine.direction.y,bisectingLine.direction.x)<=0:
			bisectingLine = Ray(PoI,PoI+bisectingLine.direction)
		else:
			bisectingLine = Ray(PoI,PoI-bisectingLine.direction)

		lIVDKey,lPoI = VDKey_intersection_at_highest_point(bisectingLine,leftVD)
		rIVDKey,rPoI = VDKey_intersection_at_highest_point(bisectingLine,rightVD)
		# print(bisectingLine)
	
	for eachKey,eachValue in lines_from_PoI_left.items():
		leftVD[eachKey]=eachValue
	for eachKey,eachValue in lines_from_PoI_right.items():
		rightVD[eachKey]=eachValue

	BisectingLines.append(((tangent.p1,tangent.p2),bisectingLine))

	VD=mergeVDUtil(leftVD,rightVD,BisectingLines)

	lowerTangent = tangent

	llen=len(leftCH)
	rlen=len(rightCH)
	ilu = leftCH.index(upperTangent.p1)
	ill = leftCH.index(lowerTangent.p1)
	iru = rightCH.index(upperTangent.p2)
	irl = rightCH.index(lowerTangent.p2)
	# print(ilu,ill,iru,irl,llen,rlen)

	# if llen>1:
	# 	while Line(leftCH[ilu],leftCH[ilu-1]).slope ==upperTangent.slope:
	# 		ilu = (ilu-1)%llen
	# 	while Line(leftCH[ill],leftCH[(ill+1)%llen]).slope ==lowerTangent.slope:	
	# 		ill = (ill+1)%llen
	
	# if rlen>1:
	# 	while Line(rightCH[iru],rightCH[(iru+1)%rlen]).slope == upperTangent.slope:	
	# 		iru = (iru+1)%rlen
	# 	while Line(rightCH[irl],rightCH[irl-1]).slope == lowerTangent.slope:		
	# 		irl = (irl-1)%rlen

	
	leftCH = leftCH[ill:ilu+1] if ilu>=ill else leftCH[ill:]+leftCH[:ilu+1]
	rightCH =rightCH[iru:irl+1] if iru<=irl else rightCH[iru:]+rightCH[:irl+1]

	return VD,leftCH+rightCH

def getVoronoiLines(points):
	points.sort(key = lambda p: [p.x,p.y])

	i=1
	while i<len(points):
		if points[i]==points[i-1]:
			print("point remove")
			points.remove(points[i]) 
		else :
			i+=1

	# print (",".join(["({0},{1})".format(p.x,p.y) for p in X]) )
	return voronoiLinesUtil(points)


if __name__=="__main__":
	# X=[-7,1,1,3,5]
	# Y=[-10,-2,3,6,-9]
	# points=[Point(a,b) for a,b in zip(X,Y)]

	npoints=random.randint(10,10)
	points=[Point(random.randint(-10,10),random.randint(-10,10)) for each in range(npoints)]
	# points=[Point(random.random(),random.random()) for each in range(npoints)]
	
	# print (",".join(["({0},{1})".format(x,y) for x,y in zip(X,Y)]) )

	VD,CH=getVoronoiLines(points)
	plot_VD_CH(VD,CH,points)
	