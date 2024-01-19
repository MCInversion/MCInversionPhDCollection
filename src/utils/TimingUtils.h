#pragma once

#include <iostream>
#include <chrono>

/**
 * \brief A simple macro for beginning of a timing session.
 * \param SESSION_NAME      chosen name (in camelCase) of the timed session.
 */
#define START_TIMING(SESSION_NAME) \
  do { \
    auto start_##SESSION_NAME = std::chrono::high_resolution_clock::now();

/**
  * \brief A simple macro for the end of a timing session.
  * \param SESSION_NAME      chosen name (in camelCase) of the timed session.
  *	\param bCoutOutput       if true, the timing output will be streamed to std::cout.
  */
#define END_TIMING(SESSION_NAME, bCoutOutput) \
    const auto end_##SESSION_NAME = std::chrono::high_resolution_clock::now(); \
    const std::chrono::duration<double> timeDiff_##SESSION_NAME = end_##SESSION_NAME - start_##SESSION_NAME; \
    if (bCoutOutput) std::cout << #SESSION_NAME ": " << timeDiff_##SESSION_NAME.count() << " s.\n"; \
  } while(0)

// Usage example:
// START_TIMING(mySession);
// /* timed code */
// END_TIMING(mySession, true);

/**
 * \brief A helper macro for timing averaged over nRuns.
 * \param SESSION_NAME         chosen name (in camelCase) of the timed session.
 * \param nRuns                the number of times the CODE_BLOCK is repeated.
 * \param CODE_BLOCK           the timed code block.
 * \param bCoutOutput          if true, the timing output will be streamed to std::cout.
 */
#define AVERAGE_TIMING(SESSION_NAME, nRuns, CODE_BLOCK, bCoutOutput) \
  do { \
    double total_duration_##SESSION_NAME = 0.0; \
    for (int run_##SESSION_NAME = 0; run_##SESSION_NAME < nRuns; ++run_##SESSION_NAME) { \
      const auto start_##SESSION_NAME = std::chrono::high_resolution_clock::now(); \
      { CODE_BLOCK } \
      const auto end_##SESSION_NAME = std::chrono::high_resolution_clock::now(); \
      total_duration_##SESSION_NAME += std::chrono::duration_cast<std::chrono::duration<double>>(end_##SESSION_NAME - start_##SESSION_NAME).count(); \
    } \
    double avg_duration_##SESSION_NAME = total_duration_##SESSION_NAME / nRuns; \
    if (bCoutOutput) std::cout << #SESSION_NAME " Average over " << nRuns << " runs: " << avg_duration_##SESSION_NAME << " s.\n"; \
  } while(0)

// Usage example:
// AVERAGE_TIMING(mySession, 10,
// {
//   /* Timed code block */
// }, true);