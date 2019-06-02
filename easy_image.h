/*
 * easy_image.h
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef EASY_IMAGE_INCLUDED
#define EASY_IMAGE_INCLUDED
#include <stdint.h>
#include <vector>
#include <iostream>
#include "Line2D.h"
#include "l_parser.h"
#include "vector3d.h"
#include "Figure.h"
#include "ZBuffer.h"
#include "Light.h"

/**
 * \brief The namespace of the EasyImage class
 */
namespace img
{
	/**
	 * \brief This class represents the color of a pixel in an img::EasyImage object
	 */
	class Color
	{
		//a safety warning to all of you: Do *NOT* rearrange the 'color components' in this class
		//easyimage expects these three fields to be the *first* fields in the class AND expects
		//them to be in *this* order
		//if you alter the arrangement, the generated BMP files will contain garbage
		public:
			/**
			 * \brief The intensity of the blue color component
			 */
			uint8_t blue;

			/**
			 * \brief The intensity of the green color component
			 */
			uint8_t green;

			/**
			 * \brief The intensity of the red color component
			 */
			uint8_t red;

			/**
			 * \brief Default Constructor
			 */
			Color();

			/**
			 * \brief Constructs a Color with the given intensities
			 *
			 * \param r	The red color component
			 * \param g	The green color component
			 * \param b	The blue color component
			 *
			 */
			Color(uint8_t r, uint8_t g, uint8_t b);

			/**
			 * Destructor
			 */
			~Color();
	};

	class Texel{
    public:
        Texel();

        Texel(int x, int y);

        int x;
	    int y;
	};

	/**
	 * \brief The exception that is thrown when an error occurs while trying to read an img::EasyImage from an input stream
	 */
	class UnsupportedFileTypeException: public std::exception
	{
		private:

			/**
			 * \brief Message explaining what went wrong
			 */
			std::string message;

		public:
			/**
			 * \brief Construct an exception with the given message
			 *
			 * \param msg	The message explaining what went wrong
			 *
			 */
			UnsupportedFileTypeException(std::string const& msg);

                        /**
                         * \brief Copy Constructor
                         *
                         * \param original	The exception to be copied into this object
                         */
			UnsupportedFileTypeException(const UnsupportedFileTypeException &original);


			/**
			 * \brief Destructor
			 */
			virtual ~UnsupportedFileTypeException() throw ();

			/**
			 * \brief Assignment operator
			 *
			 * \param original	The original exception to be assigned to this one
			 */
			UnsupportedFileTypeException& operator=(const UnsupportedFileTypeException &original);

			/**
                         * \brief Returns a description of the error hat occurred.
                         *
                         * \return A description of the error hat occurred.
                         */
			virtual const char *what() const throw ();
	};

	/**
	 * \brief This class implements a 'minor' image-library that supports basic operations such as setting and retrieving a pixel, and drawing a line.
	 */
	class EasyImage
	{
		public:
			/**
			 * \brief Default Constructor. Creates a zero-pixel image
			 */
			EasyImage();

			/**
			 * \brief Constructor: creates a new EasyImage of the specified width and height
			 *
			 * \param width		the width of the image
			 * \param height	the height of the image
			 * \param color		(optional) the background color of the image
			 */
			EasyImage(unsigned int width, unsigned int height, Color color = Color());

			/**
			 * \brief Copy Constructor
			 *
			 * \param img		the image to be copied
			 */
			EasyImage(EasyImage const& img);

			/**
			 * \brief Destructor
			 */
			virtual ~EasyImage();

			/**
			 * \brief Assignment operator. Allows an easyImage to be assigned to another easyImage
			 *
			 * \param img	The image to be assigned to this image
			 */
			EasyImage& operator=(EasyImage const& img);

			/**
			 * \brief Returns the width of the image
			 *
			 * \return the width of the image
			 */
			unsigned int get_width() const;

			/**
			 * \brief Returns the height of the image
			 * \return the height of the image
			 */
			unsigned int get_height() const;

			/**
			 * \brief Function operator. This operator returns a reference to a particular pixel of the image.
			 *
			 * \param x	the x coordinate of the pixel
			 * \param y	the y coordinate of the pixel
			 *
			 * These assertions apply:
			 * 	assert(x>=0 && x < getWidth())
			 * 	assert(y>=0 && y < getHeight())
			 */
			Color& operator()(unsigned int x, unsigned int y);

			/**
			 * \brief Function operator. This operator returns a const reference to a particular pixel of the image.
			 *
			 * \param x	the x coordinate of the pixel
			 * \param y	the y coordinate of the pixel
			 *
			 * These assertions apply:
			 * 	assert(x>=0 && x < getWidth())
			 * 	assert(y>=0 && y < getHeight())
			 */
			Color const& operator()(unsigned int x, unsigned int y) const;

			/**
			 * \brief Fills the image with a background of a specified color. Defaults to black
			 *
			 * \param color		The color to be assigned to each pixel
			 */
			void clear(Color color = Color());

			/**
			 * \brief Draws a line from pixel (x0,y0) to pixel (x1,y1) in the specified color
			 *
			 * \param x0	the x coordinate of the first pixel
			 * \param y0	the y coordinate of the first pixel
			 * \param x1	the x coordinate of the second pixel
			 * \param y1	the y coordinate of the second pixel
			 * \param color	the color of the line
			 *
			 * These assertions apply:
			 *	assert(x0 < getWidth())
			 * 	assert(y0 < getHeight())
			 * 	assert(x1 < getWidth())
			 * 	assert(y1 < getHeight())
			 */
			void draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color);

			/**
			 * \brief Draws a line from pixel (x0,y0) to pixel (x1,y1) in the specified color
			 *
			 * \param x0	the x coordinate of the first pixel
			 * \param y0	the y coordinate of the first pixel
			 * \param x1	the x coordinate of the second pixel
			 * \param y1	the y coordinate of the second pixel
			 * \param color	the color of the line
			 *
			 * These assertions apply:
			 *	assert(x0 < getWidth())
			 * 	assert(y0 < getHeight())
			 * 	assert(x1 < getWidth())
			 * 	assert(y1 < getHeight())
			 */
			void draw_zbuf_line(ZBuffer &,
								const unsigned int x0, const unsigned int y0,
								const double z0,
								const unsigned int x1, const unsigned int y1,
								const double z1,
								const Color &color);

            /**
            * \brief Draw a collection of lines
            *
            * \param Lines2D 	a list containing all the lines to draw
            * \param size       the amount of pixels on the largest side of the image
            */
            void draw2DLines(Lines2D lines2D, const int size, bool zbuffered);

            /**
            * \brief Draw a L-System
            *
            * \param l_system 	a LParser::LSystem2D containing the L-system that's going to be drawn
            * \param size   the amount of pixels on the largest side of the image
            * \param color    A vector<double> containing three numbers in [0; 1], the RGB values of the lines
            *
            * \note an l_system contains the angel in degrees, all local variables here will use the angel in radians
            */
			void drawLSystem(const LParser::LSystem2D &l_system, const int& size, std::vector<double> color);

            /**
             * \brief Generate an L-System full length string
             *
             * @param str The index string
             * @param l_system The l_system (required for the replacement rules)
             * @param size The nr of iterations for the replacement
             * @return
             */
			std::string replace_string(const std::string& str, const LParser::LSystem2D& l_system, const int& size);

			void drawFigures3D(const Figures3D &figures3D, const int size, bool zbuffered, Vector3D eyepoint);

            void draw_zbuf_triag(ZBuffer &zBuffer, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d,
                                 double dx, double dy, img::Color ambientReflection, img::Color diffuseReflection,
                                 img::Color specularReflection, double reflectionCoeff, Lights3D &lights,
                                 const Vector3D &eyepoint, bool shadow, Matrix eyePointM, Matrix invEyepointMatrix,
                                 bool texture);

            void drawZbufTriangFig(Figures3D &figures3D, const int size, Lights3D lights, const Vector3D &eyepoint,
                                   bool shadows, bool background, std::string backgroundLocation, bool texture,
                                   std::string texturefile, bool ZbufferOutput);

            void diffuseColor(img::Color &color, int xpixel, int ypixel, const double &d, const double &dx,
                              const double &dy, ZBuffer &zBuffer, const Vector3D &eyepoint, Lights3D &lights3D,
                              const Vector3D &normaalVector, const img::Color &diffuseReflection,
                              const img::Color &specularReflection, const double &reflectionCoeff, bool reflection,
                              bool shadow, Matrix eyePointM, Matrix invEyepointMatrix, bool texture);

            void applyTexture(img::EasyImage texture, double reflection);

            img::Texel generateTexel(Vector3D P);

            void setupTexture(Figures3D figures, Vector3D eyepoint, int pixelwidth, int pixelheight, double d);

		private:
			friend std::istream& operator>>(std::istream& in, EasyImage & image);
			/**
			 * \brief the width of the image
			 */
			unsigned int width;
			/**
			 * \brief the height of the image
			 */
			unsigned int height;
			/**
			 * \brief the vector containing all pixels
			 */
			std::vector<Color> bitmap;

			// Texture specific elements
			std::vector<Texel> texelmap;

			Matrix abcMatrix; // used for textures
			Vector3D p;
			int textureWidth;
			int textureHeight;

        void applyBackground(std::string backgroundLocation);
    };

	/**
	 * \brief Writes an img::EasyImage to an output stream in the BMP file format
	 *
	 * \param out		the std::ostream to write the BMP file to.
	 * \param image		the img::EasyImage to be written to the output stream
	 *
	 * \return		a reference to the output stream the image was written to
	 */
	std::ostream& operator<<(std::ostream& out, EasyImage const& image);
	/**
	 * \brief Reads an img::EasyImage from an input stream.
	 *
	 * Please note: at this point only a limited subset of BMP-file format is supported.
	 * In order to correctly read a BMP file it must:
	 * 	- Be an uncompressed bitmap
	 *	- Only contain one plane
	 * 	- Use 24bits/pixel
	 * If the BMP file-format is not supported an UnsupportedFileTypeException is thrown
	 *
	 * \param in		the input stream to read the bitmap from
	 * \param image		the EasyImage object in which the bitmap must be stored
	 *
	 * \return		a reference to the input stream from which the bitmap was read
	 */
	std::istream& operator>>(std::istream& in, EasyImage& image);

}
#endif /* EASY_IMAGE_INCLUDED */
