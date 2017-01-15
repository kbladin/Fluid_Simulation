#ifndef _APPLICATION_WINDOW_GLFW_H
#define _APPLICATION_WINDOW_GLFW_H

#include "FluidInteractionHandler.h"
#include "SGE/Controller.h"

#include <iostream>
#include <functional>

#include <gl/glew.h>
#include <gl/glfw3.h>

//! A class that handles both user interfacing and has an instance of a FluidRendererGL object.
class ApplicationWindowGLFW
{
public:
  ApplicationWindowGLFW(int width, int height);
  ~ApplicationWindowGLFW();

  void run(std::function<void(void)> f);
  void setInteractionHandler(FluidInteractionHandler* interaction_handler);
  void setController(Controller* controller);
private:
  // Functions
  bool initOpenGLContext(int width, int height);
  static void mousePosCallback(
    GLFWwindow * window,
    double x,
    double y);
  static void mouseButtonCallback(
    GLFWwindow * window,
    int button,
    int action,
    int mods);
  static void mouseScrollCallback(
    GLFWwindow * window,
    double dx,
    double dy);
  static void keyCallback(
    GLFWwindow * window,
    int key,
    int scancode,
    int action,
    int mods);
  static void windowSizeCallback(
    GLFWwindow* window,
    int width,
    int height);

  // Data
  GLFWwindow* _window;
  static FluidInteractionHandler* _interaction_handler;
  static Controller* _controller;

  float _delay_counter;
  int   _frame_counter;
  double _time;
};

#endif
