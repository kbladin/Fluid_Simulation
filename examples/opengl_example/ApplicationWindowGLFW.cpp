#include "ApplicationWindowGLFW.h"

FluidInteractionHandler* ApplicationWindowGLFW::_interaction_handler;

ApplicationWindowGLFW::ApplicationWindowGLFW(int width, int height)
{
  _frame_counter = 0;
  _delay_counter = 0;
  _time = glfwGetTime();
  // First an OpenGL context needs to be created
  if (!initOpenGLContext(width, height))
  {
    std::cout << "ERROR : Failed to initialize OpenGL" << std::endl;
  }
  glfwSwapInterval(1);
  // Set callback functions
  glfwSetCursorPosCallback(_window, mousePosCallback);
  glfwSetMouseButtonCallback(_window, mouseButtonCallback);
  glfwSetScrollCallback(_window, mouseScrollCallback);
  glfwSetKeyCallback(_window, keyCallback);
  glfwSetWindowSizeCallback(_window, windowSizeCallback);
}

ApplicationWindowGLFW::~ApplicationWindowGLFW()
{
  glfwTerminate();
}

bool ApplicationWindowGLFW::initOpenGLContext(int width, int height)
{
  // Initialize glfw
  if (!glfwInit())
    return false;
  // Modern OpenGL
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  // Create a windowed mode window and its OpenGL context
  _window = glfwCreateWindow(width, height, "Fluid Simulation", NULL, NULL);
  if (!_window)
  {
    glfwTerminate();
    return false;
  }
  // Make the window's context current
  glfwMakeContextCurrent(_window);
  printf("%s\n", glGetString(GL_VERSION));
  return true;
}

//! Starts the main loop
void ApplicationWindowGLFW::run(std::function<void(void)> f)
{
  while (!glfwWindowShouldClose(_window))
  {
    f();

    _frame_counter ++;
    double time_since_last = glfwGetTime() - _time;
    _delay_counter += time_since_last;

    if (_delay_counter >= 1.0) {
      std::stringstream title;
      title << "Fluid Simulation. " << _frame_counter << " FPS";
      glfwSetWindowTitle(_window, title.str().c_str());
      _frame_counter = 0;
      _delay_counter = 0;
    }


    glfwSwapBuffers(_window);
    glfwPollEvents();
    _time = glfwGetTime();
  }
}

void ApplicationWindowGLFW::setInteractionHandler(FluidInteractionHandler* interaction_handler)
{
  _interaction_handler = interaction_handler;
}

void ApplicationWindowGLFW::mousePosCallback(
  GLFWwindow * window,
  double x,
  double y)
{
  if (_interaction_handler)
  {
    int size_x, size_y;
    glfwGetWindowSize(window, &size_x, &size_y);
    _interaction_handler->mousePosCallback(
      (x / size_x - 0.5) * 2,
      (1 - (y / size_y) - 0.5) * 2);
  }
}

void ApplicationWindowGLFW::mouseButtonCallback(
  GLFWwindow * window,
  int button,
  int action,
  int mods)
{
  if (_interaction_handler)
  {
    _interaction_handler->mouseButtonCallback(button, action, mods);
  }
}

void ApplicationWindowGLFW::mouseScrollCallback(
  GLFWwindow * window,
  double dx,
  double dy)
{
  if (_interaction_handler)
  {
    _interaction_handler->mouseScrollCallback(dx, dy);
  }
}

void ApplicationWindowGLFW::keyCallback(
  GLFWwindow * window,
  int key,
  int scancode,
  int action,
  int mods)
{
  if (_interaction_handler)
  {
    _interaction_handler->keyCallback(key, scancode, action, mods);
  }
}

void ApplicationWindowGLFW::windowSizeCallback(
  GLFWwindow* window,
  int width,
  int height)
{
  if (_interaction_handler)
  {
    _interaction_handler->windowSizeCallback(width, height);
  }
}
