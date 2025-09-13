#include <catch2/catch_all.hpp>
#include <chrono>
#include <thread>
#include <iostream>
#include "ControlCHandler.h"

using namespace RDKit;
using namespace std::chrono_literals;

int main(int argc, char *argv[]) {
  ControlCHandler::reset();
  for (;;) {
    if (ControlCHandler::getGotSignal()) {
      std::cerr << "CTRL+C was pressed" << std::endl;
      break;
    }
    std::this_thread::sleep_for(10ms);
  }
  return 0;
}
