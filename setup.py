from distutils.core import setup

setup(name="pyurdme",
      version="1.0.0",
      author="Andreas Hellander, Brian Drawert",
      author_email="andreas.hellander@gmail.com",
      packages=['pyurdme','examples','examples/mincde','examples/cylinder_demo'],
      package_data={'examples/mincde':['mesh/*']})
