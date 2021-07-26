#ifndef FRAMESTEPPER_H
#define FRAMESTEPPER_H

class FrameStepper {

   float frame_length; //length of a frame in seconds
   int frame_count; //always the frame currently being processed
   float current_time; //current time within a frame
   int step_count; //current step within the frame
   int halt_frame; //number of frames for the simulation

public:
   
   FrameStepper (float frame_len, int frame_limit);
   
   //adjust the timestep to land on a frame time, or to use more evenly spaced steps if close to a frame time
   float get_step_length(float max_step);

   //we're done when current time is very close or past the frame_length
   bool done_frame();

   //sim is complete when we hit our target frame_count
   bool done_simulation();

   void advance_step(float step_length);
   void advance_frame();

   int get_step_count();
   int get_frame();
   float get_time();
};

#endif
