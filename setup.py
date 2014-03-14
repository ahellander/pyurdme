from distutils.core import setup

setup(name="pyurdme",
      version="0.1.0",
      author="Andreas Hellander, Brian Drawert",
      author_email="andreas.hellander@gmail.com",
      packages=['pyurdme'],
      package_data={'pyurdme':['data/*.c','data/three.js_templates/*','urdme/AUTHORS','urdme/LICENCE','urdme/bin/*','urdme/build/*','urdme/include/*','urdme/src/*.c','urdme/src/nsm/*']})