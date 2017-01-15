#include "ApplicationWindowGLFW.h"

FluidInteractionHandler* ApplicationWindowGLFW::_interaction_handler;
Controller* ApplicationWindowGLFW::_controller;

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
  windowSizeCallback(_window, 0, 0); // Just to refresh
}

void ApplicationWindowGLFW::setController(Controller* controller)
{
  _controller = controller;
}

void ApplicationWindowGLFW::mousePosCallback(
  GLFWwindow * window,
  double x,
  double y)
{
  int size_x, size_y;
  glfwGetWindowSize(window, &size_x, &size_y);
  // Convert to NDC
  x = (x / size_x - 0.5) * 2;
  y = (1 - (y / size_y) - 0.5) * 2;
  if (_interaction_handler)
  {
    _interaction_handler->mousePosCallback(x, y);
  }
  if (_controller)
  {
    _controller->mousePosCallback(x, y);
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
  if (_controller)
  {
    _controller->mouseButtonCallback(
      static_cast<MouseButton>(button), static_cast<KeyAction>(action));
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
  if (_controller)
  {
    _controller->mouseScrollCallback(dx, dy);
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
  if (_controller)
  {
    _controller->keyCallback(
      static_cast<Key>(key), static_cast<KeyAction>(action));
  }
}

void ApplicationWindowGLFW::windowSizeCallback(
  GLFWwindow* window,
  int width,
  int height)
{
  int frame_buffer_size_x, frame_buffer_size_y;
  glfwGetFramebufferSize(window, &frame_buffer_size_x, &frame_buffer_size_y);
  if (_interaction_handler)
  {
    _interaction_handler->windowSizeCallback(frame_buffer_size_x, frame_buffer_size_y);
  }
  if (_controller)
  {
    _controller->windowSizeCallback(frame_buffer_size_x, frame_buffer_size_y);
  }
}
