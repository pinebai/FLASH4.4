class TimeDep:
    def __init__(self, fig, num_frames, init_frame, paint):

        if num_frames() == 0:
            raise ValueError("number of frames must be greater than zero")

        self.current_frame = init_frame
        self.num_frames = num_frames
        self.fig = fig
        self.paint = paint

        cid = self.fig.canvas.mpl_connect('key_press_event', self.onkey)

        self.paint(init_frame, self.fig)

    def onkey(self, event):
        if event.key == 'right' and self.current_frame != self.num_frames()-1:
            self.current_frame += 1
            self.paint(self.current_frame, self.fig)

        if event.key == 'left'  and self.current_frame != 0:
            self.current_frame -= 1
            self.paint(self.current_frame, self.fig)

        if event.key == 'pageup' and self.current_frame != self.num_frames()-1:

            self.current_frame += 10
            if self.current_frame > self.num_frames()-1: 
                self.current_frame = self.num_frames()-1

            self.paint(self.current_frame, self.fig)


        if event.key == 'pagedown' and self.current_frame != 0:

            self.current_frame -= 10
            if self.current_frame < 0:
                self.current_frame = 0

            self.paint(self.current_frame, self.fig)
