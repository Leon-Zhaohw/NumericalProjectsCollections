#include "framestepper.h"

FrameStepper::FrameStepper (float frame_len, int frame_limit):
frame_length(frame_len), 
frame_count(1), 
step_count(1),
halt_frame(frame_limit)
{}

//adjust the timestep to land on a frame time, or to use more evenly spaced steps if close to a frame time
float FrameStepper::get_step_length(float max_step) {
   if(current_time + max_step > frame_length)
      max_step = frame_length - current_time;
   else if(current_time + 1.5f*max_step >= frame_length)
      max_step = 0.5f*(frame_length-current_time);

   return max_step;     
}

//we're done when current time is very close or past the frame_length
bool FrameStepper::done_frame() {
   return current_time >= 0.99*frame_length;
}

//sim is complete when we hit our target frame_count
bool FrameStepper::done_simulation() {
   return frame_count == halt_frame;
}

void FrameStepper::advance_step(float step_length) {
   current_time += step_length;
   ++step_count;
}

void FrameStepper::advance_frame() {
   current_time = 0;
   step_count = 1;
   ++frame_count;
}

int FrameStepper::get_step_count() {
   return step_count;
}

int FrameStepper::get_frame() {
   return frame_count;
}

float FrameStepper::get_time() {
   return (frame_count-1)*frame_length + current_time;

}