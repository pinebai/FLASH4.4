 #!usr/bin/python
from Tkinter import *
from math import *
import thread, time, string, struct, array, os



#This routine will plot particles
#usage > python particles.py test_DumpParticles_0000 numParticles

print sys.argv




class gui:
  width = 700
  height = 700
  xmax = 10.
  xmin = 0.
  ymax = 10.
  ymin = 0.
  part_sz = 3.0  
  sleep = 0.01
  trace = 0
  

  def screen_coords(self, x, y):
    pos = []
    pos.append(((x-gui.xmin)/(gui.xmax-gui.xmin))*gui.width)
    pos.append(((y-gui.ymin)/(gui.ymax-gui.ymin))*gui.height)
    return pos



  def run(self):
    #input, num_particles = gui.open_file(self)
    particles_data, num_particles = gui.read_binary2(self)
    gui.plot(self,particles_data, num_particles)
    #gui.close_file(self, input)

  def open_file(self):
    filename = sys.argv[1]
    print filename
    num_particles = int(sys.argv[2])
    print num_particles
    input = open(filename, 'r')
    #gui.read(self, input, num_particles)
    #gui.assemble(self, scaf_input, target_input, 1, 1, 1)
    return input, num_particles


  def close_file(self, input):
    input.close()
    

  def read_binary(self):
    filename = sys.argv[1]
    print filename
    num_particles = int(sys.argv[2])
    print num_particles
    input = open(filename, 'rb')
    
    particle_data = [] #holds all particles, timesteps

    particleArr = array.array("f") #f for 'f'loating point
    numItems = os.stat(filename)[6]/particleArr.itemsize
    particleArr.fromfile(open(filename, "rb"), numItems)
    particle_data.append(particleArr)
                      
#    print "particle_data = ", particle_data                  
  

  def read_binary2(self):
    filename = sys.argv[1]
    print filename
    num_particles = int(sys.argv[2])
    print num_particles

    
    data = open(filename).read()
    start,stop = 0,struct.calcsize('d')
    print "size = ", stop

    datalen = len(data)
    print "datalen = ", datalen

    size = stop

    particle_struct = []

    while stop <= datalen:
      particle_data = []
      for i in range(num_particles):
        single_particle = []
        for j in range(4):  #4 because storing tag, blk, xpos and ypos
          part = struct.unpack('d', data[start:stop])
          start = stop
          stop = stop + size
          single_particle.append(part[0])

        particle_data.append(single_particle)  
      particle_struct.append(particle_data)  

    #print "particle_struct = ", particle_struct  

    return particle_struct, num_particles
    

  def read(self, input, num_particles):
    particle_struct = [] #holds all particles, timesteps

    not_end_of_file = 1
    
    while not_end_of_file:
      single_particle= []
      num_particle_data = []

      for j in range(num_particles):
        line = input.readline()
        if not line:
          not_end_of_file = 0
          break
        particle = [ x.strip() for x in line.split() ]
        num_particle_data.append(particle)

      particle_struct.append(num_particle_data)  

    for i in range(10):
      print 'step = ', i
      print 'particle 1'
      print particle_struct[i][0]
      print 'particle 2'
      print particle_struct[i][1]

    return particle_struct
      


  def plot(self, particles_data, num_particles):
    print len(particles_data)


    #create the initial positions of the particles.  
    pos = []      
    part = []
    for j in range(num_particles):
      
      orig_pos = gui.screen_coords(self,float(particles_data[0][j][2]), float(particles_data[0][j][3]))
      pos.append(orig_pos)
#      print "orig xpos, orig ypos ", orig_pos[0], orig_pos[1]
    
      color = 'black'
      ocolor = 'black'
      
      particle = (c.create_oval(pos[j][0] - self.part_sz, pos[j][1] - self.part_sz, pos[j][0] + 
                                self.part_sz, pos[j][1] + self.part_sz, fill=color, outline=ocolor))
      
      part.append(particle)
      



    #now either move the particles or leave a trace (meaning draw another circle)    
    for i in range(len(particles_data) -2):
    
      new_pos = []      
      for j in range(num_particles):


#        print "raw x, raw y ", particles_data[i+1][j][2], particles_data[i+1][j][2]
        npos = gui.screen_coords(self,float(particles_data[i+1][j][2]), float(particles_data[i+1][j][3]))
        new_pos.append(npos)
        
        if gui.trace:

          
          part.append(c.create_oval(new_pos[j][0] - self.part_sz, new_pos[j][1] - self.part_sz, new_pos[j][0] + 
                                    self.part_sz, new_pos[j][1] + self.part_sz, fill=color, outline=ocolor))


        else:  #don't do tracing, just move the original particle
      
                    
          #print 'step = ', str(i)
          #print 'particle = ', str(part[j])
          #print 'old pos = ', str(pos)
          #print 'new pos = ', str(new_pos)
          
          c.move(part[j], new_pos[j][0] - pos[j][0], new_pos[j][1] - pos[j][1])
#          print "xpos, ypos ", new_pos[j][0], new_pos[j][1]
          pos[j] = new_pos[j]
          
        
      time.sleep(gui.sleep)
      root.update()






#ge = gui()
#ge.read_binary2()




ge = gui() 





root = Tk()
root.title('ViewParticles')
root.resizable(width=NO, height=NO)
menu_win = Canvas(root, width=ge.width, height=50, background='white')
c = Canvas(root, width=ge.width, height=ge.height, background='white')
c.config(scrollregion=(0,0, 200*ge.width, 30*ge.height))
c.config(highlightthickness=0)
vbar = Scrollbar(root)
vbar.config(command=c.yview)
c.config(yscrollcommand=vbar.set)
vbar.pack(side=RIGHT, fill=Y)

hbar = Scrollbar(root, orient='horizontal')
hbar.config(command=c.xview)
c.config(xscrollcommand=hbar.set)
hbar.pack(side=BOTTOM, fill=X)

run_but = Button(menu_win, text= 'run', command=ge.run)
run_but.grid(row=0, column=0)
exit_but = Button(menu_win, text= 'exit', command=root.quit)
exit_but.grid(row=0, column=5)
menu_win.pack()
c.pack()
c.focus()
root.mainloop()
