#!/usr/local/bin/visit_script


from sys import exit
from os import access, R_OK

from math import pi

path = "/home/ahaque/THREE/trunk/SB2S"
data_subdir = "dumps"
output_subdir = "dumps"
sim_name = "sb"
data_name = sim_name + "_hdf5_chk_"

data_file = "data"
data_suffix = "curve"

temp_path = path + "/" + data_subdir
temp_file = temp_path + "lineout.temp"

file_index = int(0)
index_digits = int(4)

var = "vely"

start_index = int(0)
end_index = int(100)

num_samples = 10000

radius = float(1.0)

## Python 2.1.2
True = int(1)
False = int(0)
## End Python 2.1.2

from sys import stdout, stderr
from os import system


num_angles = 3
dtheta = pi / 4.0


# OpenComputeEngine("fornax", ("-np", "2"))


class Database :
    def __init__(self, radius, num_samples, variable="ch  ",
                 temp_lineout_name="/tmp/lineout.temp") :

        self.__radius = radius

        self.__num_samples = num_samples

        self.__variable = variable

        self.__temp_lineout_name = temp_lineout_name

        self.__dbase_open = False
        self.__database_filename = None


    def open_database(self, filename) :
        if (self.__dbase_open) :
            self.close_database()
        
        OpenDatabase(filename)

        AddPlot("Pseudocolor", self.__variable)
        DrawPlots()

        self.__dbase_open = True
        self.__database_filename = filename


    def close_database(self) :
        if (self.__dbase_open) :
            DeleteAllPlots()
            # DeleteWindow()
            CloseDatabase(self.__database_filename)

            self.__database_filename = None

        self.__dbase_open = False


    def get_time(self) :
        if (not self.__dbase_open) :
            stderr.write("No database open\n")
            raise str("No database open")

        Query("Time")

        return GetQueryOutputValue()


    def get_line(self, angle, next_window_id=2) :
        from math import sin, cos
        
        if (not self.__dbase_open) :
            stderr.write("No database open\n")
            raise str("No database open")

        # Produce the lineout plot


        Lineout((0.0, 0.0), ((self.__radius * cos(angle)),(self.__radius * sin(angle))),
                self.__num_samples)


        # Save the data

        SetActiveWindow(next_window_id)

        sw_attr = GetSaveWindowAttributes()
        sw_attr.fileName = self.__temp_lineout_name
        sw_attr.format = sw_attr.CURVE

        SetSaveWindowAttributes(sw_attr)

        filename = SaveWindow()

        DeleteWindow()

        # Read the data

        radii = [ ]
        values = [ ]

        f = open(filename, "r")

        file_lines = f.readlines()

        for line in file_lines :
            words = line.split()

            num_words = len(words)
            
            if (num_words == 2) :
                if (words[0][0] != '#') :
                    radii.append(float(words[0]))
                    values.append(float(words[1]))

        f.close()

        system("rm -f " + filename)

        return radii, values


def build_database_name(index) :
    return path + "/" + data_subdir + "/" + data_name \
           + ("%(idx)04d" % { "idx" : index })


def main() :

    dbase = Database(radius, num_samples, var,
                     temp_lineout_name=temp_file)

    for current_file_index in range(start_index, (end_index+1)) :
        
        # Open the database

        current_file = build_database_name(current_file_index)
    
        dbase.open_database(current_file)

        stdout.write("Processing database [ " + current_file + " ]\n")

        time = dbase.get_time()

        print time

        for angle_index in range(num_angles) :

            angle = float(angle_index) * dtheta

            # Open the output file

            output_filename = path + "/" + output_subdir + "/" + sim_name \
                              + "_" + var + "_" + str(current_file_index) + "_" + str(angle_index) + ".out"

            f = open(output_filename, "w+")

            stdout.write("Opened file [ " + output_filename + " ]\n")

            radii, values = dbase.get_line(angle)

            num_points = len(radii)

            f.write("# "+str(time)+"\n")
            f.write("# Angle " + str(angle) + "\n")
            f.write("# Flash file index " + str(current_file_index) + "\n")

            for i in range(num_points) :
                f.write(str(radii[i]) + " " + str(values[i]) + "\n")


if (__name__ == "__main__") :
    
    main()
    exit()
    
