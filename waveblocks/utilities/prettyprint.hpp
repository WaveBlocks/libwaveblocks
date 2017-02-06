#pragma once

namespace waveblocks {
namespace utilities {

/**
 * \file
 * \namespace prettyprint
 * \brief simple helper functions for prettier command line output
 */
namespace prettyprint {

	static unsigned WIDTH = 60; ///< width of the command line output

	/**
	 * \brief prints a string-value pair
	 * Prints a string-value pair that is aligned left and right of the available space respectively
	 * \param s string (aligned on the left side of the printed output)
	 * \param v value (aligned on the right side of the printed output)
	 * \param s0 optional argument that is written to std::cout before printing the pair.
	 *        The default value is a line break "\n".
	 *        If the previous line should be overwritten, set s0="\r"
	 */
	template <typename T, typename S>
	inline void pair(const T key, const S val, const std::string&& s0="\n") {
		std::streamsize default_precision = std::cout.precision();
		std::cout << s0 << "\t"
			<< std::setw(WIDTH/2) << std::left << key
			<< std::setw(WIDTH/2) << std::right << std::fixed << std::setprecision(4) << val
			<< std::scientific << std::setprecision(default_precision) << std::flush;
	}
	
	/**
	 * \brief prints a string as a title
	 * Prints a centered, decorated string to std::cout
	 * \param s string to be used for the title
	 * \param c optional delimiting character used for decorating the title
	 *        The default value is a simple dash '-'.
	 */
	inline void title(const std::string&& s, const char c='-') {
		int str_length = s.length();
		int w1 = (WIDTH-str_length-2)/2;
		int w2 = (WIDTH-str_length-2) - w1;
		std::cout
			<< "\n\t" << std::string(WIDTH,c)
			<< "\n\t" << c << std::string(w1,' ') << s << std::string(w2,' ') << c
			<< "\n\t" << std::string(WIDTH,c);
	}

	/**
	 * \brief prints a horizontal separator
	 * \param c optional delimiting character to be used for the separator.
	 *        The default value is a simple dash '-'.
	 */
	inline void separator(const char c='-') {
		std::cout << "\n\t" << std::string(WIDTH,c) << std::flush;
	}

} // namespace prettyprint
} // namespace utilities
} // namespace waveblocks
