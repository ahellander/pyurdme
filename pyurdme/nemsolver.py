""" NSM solver. """
import pyurdme
import os

class NEMSolver(pyurdme.URDMESolver):
    """ NEM solver class. """

    NAME = 'nem'

    def createPropensityFile(self, file_name=None):
        """ Generate a C propensity file on the new experimental format. """

        template = open(os.path.abspath(os.path.dirname(__file__)) + '/data/propensity_file_new_template.c', 'r')
        propfile = open(file_name, "w")
        propfilestr = template.read()

        propfilestr = propfilestr.replace("__NUMBER_OF_REACTIONS__", str(self.model.getNumReactions()))
        propfilestr = propfilestr.replace("__NUMBER_OF_SPECIES__", str(len(self.model.listOfSpecies)))


        speciesdef = ""
        for i, sname in enumerate(self.model.listOfSpecies):
            S = self.model.listOfSpecies[sname]
            speciesdef += "species *" + sname + ";\n\t"
            speciesdef += sname + "= (species *)malloc(sizeof(species));\n\t"
            speciesdef += sname + "->gamma = " + str(S.diffusion_constant) + ";\n\t"
            speciesdef += sname + "->sigma = " + str(S.reaction_radius) + ";\n\t"
            speciesdef += "ptr["+str(i)+"] = " + sname + ";\n\n\t"


        propfilestr = propfilestr.replace("__DEFINE_SPECIES__", speciesdef)


        # Make sure all paramters are evaluated to scalars before we write them to the file.
        self.model.resolveParameters()

        reacstr = ""

        for j, sname in enumerate(self.model.listOfSpecies):
            reacstr += "int " + sname + "=" + str(j) + ";\n\t"

        reacstr += "\n\t"

        for i, rname in enumerate(self.model.listOfReactions):
            R = self.model.listOfReactions[rname]

            reacstr += "reaction *" + rname + ";\n\t"
            reacstr += rname + "=(reaction *)malloc(sizeof(reaction));\n\t"
            reacstr += rname + "->order=" + str(len(R.reactants)) + ";\n\t"
            reacstr += rname + "->nr_reactants=" + str(len(R.reactants)) + ";\n\t"
            reacstr += rname + "->nr_products=" + str(len(R.products)) + ";\n\t"


            reacstr += rname+"->reactants=(int *)malloc(" + rname + "->nr_reactants*sizeof(int));\n\t"
            for j, reactant in enumerate(R.reactants):
                reacstr += rname + "->reactants[" + str(j) + "]=" + str(reactant)+";\n\t"

            reacstr += "\n\t" + rname + "->products=(int *)malloc(" + rname + "->nr_products*sizeof(int));\n\t"
            for j, product in enumerate(R.products):
                reacstr += rname + "->products[" + str(j) + "]=" + str(product) + ";\n\t"

            reacstr += "\n\t" + rname + "->nr=(int *)calloc(" + str(len(self.model.listOfSpecies)) + ",sizeof(int));\n\t"

            for j, reactant in enumerate(R.reactants):
                reacstr += rname + "->nr[" + reactant + "]=-" + str(R.reactants[reactant]) + ";\n\t"

            for j, product in enumerate(R.products):
                reacstr += rname + "->nr[" + product + "]=" + str(R.products[product]) + ";\n\t"

            reacstr += rname + "->k=" + str(R.marate.value) + ";\n\t"

            reacstr += "\n\tptr[" + str(i) + "] = " + rname + ";\n\n\t"


        propfilestr = propfilestr.replace("__DEFINE_REACTIONS__", reacstr)

        propfile.write(propfilestr)
        propfile.close()
