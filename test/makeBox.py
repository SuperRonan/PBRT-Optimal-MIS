

size = 3

useSpehre = True

dark = 0.1
material = 'MakeNamedMaterial "TheDarkestMaterialEver" "string type" [ "matte" ] "rgb Kd" [ %f %f %f ]' % (dark, dark, dark)


res = 'NamedMaterial "TheDarkestMaterialEver"'

if(useSpehre):
    res = res + '''
Shape "sphere"
"float radius" [%f] ''' % (size)
else:
    res = res + '''
    Scale %f %f %f 
    TransformBegin
# floor
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ -1 -1 -1   -1 -1 1   1 -1 1   1 -1 -1 ] "normal N" [ 0 1 0   0 1 0   0 1 0   0 1 0 ] "float uv" [ 0 0 1 0 1 1 0 1 ] 
	
# ceiling
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ 1 1 1   -1 1 1   -1 1 -1   1 1 -1 ] "normal N" [ 0 -1 0   0 -1 0   0 -1 0   0 -1 0 ] "float uv" [ 0 0 1 0 1 1 0 1 ]

# wall front of the camera
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ -1 -1 -1   -1 1 -1   1 1 -1   1 -1 -1 ] "normal N" [ 0 0 -1   0 0 -1   0 0 -1   0 0 -1 ] "float uv" [ 0 0 1 0 1 1 0 1 ] 
	
# wall behind the camera
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ -1 -1 1   -1 1 1   1 1 1   1 -1 1 ] "normal N" [ 0 0 1   0 0 1   0 0 1   0 0 1 ] "float uv" [ 0 0 1 0 1 1 0 1 ] 
	
# right wall
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ 1 -1 -1   1 1 -1   1 1 1   1 -1 1 ] "normal N" [ 1 0 0   1 0 0   1 0 0   1 0 0 ] "float uv" [ 0 0 1 0 1 1 0 1 ] 
	
# left wall
	Shape "trianglemesh" "integer indices" [ 0 1 2 0 2 3 ] "point P" [ -1 -1 1   -1 1 1   -1 1 -1   -1 -1 -1 ] "normal N" [ -1 0 0   -1 0 0   -1 0 0   -1 0 0 ] "float uv" [ 0 0 1 0 1 1 0 1 ] 
    TransformEnd
    ''' % (size, size, size)

#res = "AttributeBegin\n" + res + "\nAttributeEnd"

print(material)
print(res)