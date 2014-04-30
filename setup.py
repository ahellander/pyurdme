from setuptools import setup

setup(name="pyurdme",
      version="1.0.0-c1",
      packages=['pyurdme'],
      
      #include_package_data = True,
      package_data={'pyurdme':['data/*.c','data/three.js_templates/js/*','data/three.js_templates/*.html','urdme/AUTHORS','urdme/LICENCE','urdme/bin/*','urdme/build/*','urdme/include/*','urdme/src/*.c','urdme/src/nsm/*']},
      
      install_requires = ["numpy","matplotlib","scipy","h5py","Jinja2"],
      
      author="Andreas Hellander, Brian Drawert",
      author_email="andreas.hellander@gmail.com",
      license = "GPL",
      keywords = "pyurdme, urdme, spatial stochastic simulation, rdme",
      #url = ""
      
      )
      
