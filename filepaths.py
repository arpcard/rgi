import sys
import os

'''
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
'''
def determine_path():
    try:
        root = __file__
        if os.path.islink (root):
            root = os.path.realpath (root)
        return os.path.dirname (os.path.abspath (root))
    except:
        print "I'm sorry, but something is wrong."
        print "There is no __file__ variable. Please contact the author."
        sys.exit ()

#print determine_path()