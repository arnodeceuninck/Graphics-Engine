#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "Figure.h"
#include "TransMatrix.h"
#include "3DPlatonic.h"
#include "LSystem3D.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>


/*
 * Graphics Engine - Arno Deceuninck
 *
 * Artefacten?
 *      - Overal double gebruiken: check
 *      - Enkel vlak voor het tekenen afronden: check
 *      - roundToInt (roundInt) functie gebruiken: check
 *      => Waarschijnlijk door klein verschil in berekeningen tov referentie-implementatie
 *
 * Dubbele lijnen?
 *      // TODO: zorgen dat eens er een lijn van punt A naar B getekend is, dit niet nog eens getekend word
 *
 * Correct swappen z-waarden bij z buffering?
 *      // TODO
 *
 *  // TODO: Stochastische L systemen
 *  // TODO: Menger sponge bij schaduwen vreemd probleem
 *  // TODO: punt fout bij bucky ball (zichtbaar bij spheres and cylinders)
 *  // TODO: L systemen bij spheres and cylinders op stukken die rechtoot gaan
 *  // TODO: Fix zbuffer als je er recht op kijkt?
 */

img::EasyImage coloredRectangle() {
    img::EasyImage image(256, 256);
    for (unsigned int i = 0; i < 256; i++) {
        for (unsigned int j = 0; j < 256; j++) {
            image(i, j).red = i;
            image(i, j).green = j;
            image(i, j).blue = (i + j) % 256;
        }
    }
    return image;
}

img::EasyImage
l_system_2D(std::vector<double> backgroundcolor, std::vector<double> color, std::string inputfile, int size) {
    LParser::LSystem2D l_system;

    std::ifstream input_stream(inputfile);
    input_stream >> l_system;
    input_stream.close();

    img::EasyImage image(size, size,
                         img::Color(backgroundcolor[0] * 255, backgroundcolor[1] * 255, backgroundcolor[2] * 255));
    image.drawLSystem(l_system, size, color);

    return image;
}

Figure l_system_3D(std::vector<double> backgroundcolor, std::vector<double> color, std::string inputfile) {
    LParser::LSystem3D l_system;

    std::ifstream input_stream(inputfile);
    input_stream >> l_system;
    input_stream.close();

    Figure image = draw3DLSystem(l_system, color);

    return image;
}

img::EasyImage generate_image(const ini::Configuration &configuration) {
    // Bepaal het type van de afbeelding
    std::string type = configuration["General"]["type"].as_string_or_die();

    // parameters nodig voor alle afbeeldingen
    int size = configuration["General"]["size"].as_int_or_default(512);
    std::vector<double> backgroundcolor = configuration["General"]["backgroundcolor"].as_double_tuple_or_default(
            {0, 0, 0});
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_default(0);
    std::vector<double> eye_vector = configuration["General"]["eye"].as_double_tuple_or_default({0, 0, 0});

    bool shadowEnabeled = configuration["General"]["shadowEnabled"].as_bool_or_default(false);

    bool zBufferOutput = configuration["General"]["shadowEnabled"].as_bool_or_default(false);
    // Zoek de bijhorende functie
    // Haal alle parameters voor bij die functie en return dan ook die afbeelding
    if (type == "2DLSystem") {
        std::string inputfile = configuration["2DLSystem"]["inputfile"].as_string_or_die();
        std::vector<double> color = configuration["2DLSystem"]["color"].as_double_tuple_or_default({1, 1, 1});
        return l_system_2D(backgroundcolor, color, inputfile, size);
    } else if (type == "Wireframe" or type == "ZBufferedWireframe" or type == "ZBuffering" or
               type == "LightedZBuffering") {
        bool light = type == "LightedZBuffering";
        Figures3D figures3D = {};
        for (int i = 0; i < nrFigures; ++i) {

            bool fractal = false;
            bool thickFigure = false;

            std::string figurename = "Figure" + std::to_string(i);

            std::string figure_type = configuration[figurename]["type"].as_string_or_die();

            Figure figure = Figure();

            if (figure_type == "LineDrawing" or figure_type == "ThickLineDrawing") {
                if (figure_type == "ThickLineDrawing") {
                    thickFigure = true;
                }

                int nrPoints = configuration[figurename]["nrPoints"].as_int_or_default(0);
                int nrLines = configuration[figurename]["nrLines"].as_int_or_default(0);

                for (int j = 0; j < nrPoints; ++j) {
                    std::string pointname = "point" + std::to_string(j);
                    std::vector<double> point_vec = configuration[figurename][pointname].as_double_tuple_or_die();
                    figure.points.emplace_back(Vector3D::point(point_vec[0], point_vec[1], point_vec[2]));
                }

                for (int k = 0; k < nrLines; ++k) {
                    std::string linename = "line" + std::to_string(k);
                    std::vector<int> line_vec = configuration[figurename][linename].as_int_tuple_or_die();
                    figure.faces.emplace_back(Face(line_vec));
                }

            } else if (figure_type == "ThickCube") {
                thickFigure = true;
                figure = createCube();
            } else if (figure_type == "ThickTetrahedron") {
                thickFigure = true;
                figure = createTetrahedron();
            } else if (figure_type == "ThickOctahedron") {
                thickFigure = true;
                figure = createOctahedron();
            } else if (figure_type == "ThickIcosahedron") {
                thickFigure = true;
                figure = createIcosahedron();
            } else if (figure_type == "ThickDodecahedron") {
                thickFigure = true;
                figure = createDodecahedron();
            } else if (figure_type == "ThickBuckyBall") {
                thickFigure = true;
                figure = createBuckeyBall();
            } else if (figure_type == "FractalCube") {
                fractal = true;
                figure = createCube();
            } else if (figure_type == "FractalTetrahedron") {
                fractal = true;
                figure = createTetrahedron();
            } else if (figure_type == "FractalOctahedron") {
                fractal = true;
                figure = createOctahedron();
            } else if (figure_type == "FractalIcosahedron") {
                fractal = true;
                figure = createIcosahedron();
            } else if (figure_type == "FractalDodecahedron") {
                fractal = true;
                figure = createDodecahedron();
            } else if (figure_type == "FractalBuckyBall") {
                fractal = true;
                figure = createBuckeyBall();
            } else if (figure_type == "Cube") {
                figure = createCube();
            } else if (figure_type == "Tetrahedron") {
                figure = createTetrahedron();
            } else if (figure_type == "Octahedron") {
                figure = createOctahedron();
            } else if (figure_type == "Icosahedron") {
                figure = createIcosahedron();
            } else if (figure_type == "Dodecahedron") {
                figure = createDodecahedron();
            } else if (figure_type == "Cylinder") {
                double height = configuration[figurename]["height"].as_double_or_default(5);
                int n = configuration[figurename]["n"].as_int_or_default(20);
                figure = createCylinder(n, height, false);
            } else if (figure_type == "Cone") {
                double height = configuration[figurename]["height"].as_double_or_default(5);
                int n = configuration[figurename]["n"].as_int_or_default(20);
                figure = createCone(n, height);
            } else if (figure_type == "Sphere") {
                int n = configuration[figurename]["n"].as_int_or_default(20);
                figure = createSphere(n);
            } else if (figure_type == "Torus") {
                double r = configuration[figurename]["r"].as_double_or_default(20);
                double R = configuration[figurename]["R"].as_double_or_default(20);
                int m = configuration[figurename]["m"].as_int_or_default(20);
                int n = configuration[figurename]["n"].as_int_or_default(20);
                figure = createTorus(r, R, n, m);
            } else if (figure_type == "3DLSystem" or figure_type == "Thick3DLSystem") {
                thickFigure = figure_type == "Thick3DLSystem";
                std::string inputfile = configuration[figurename]["inputfile"].as_string_or_die();
                std::vector<double> color = configuration[figurename]["color"].as_double_tuple_or_default({1, 1, 1});
                figure = l_system_3D(backgroundcolor, color, inputfile);
            } else if (figure_type == "BuckyBall") {
                figure = createBuckeyBall();
            } else if (figure_type == "MengerSponge") {
//                std::cout << "WIP, please come back later" <<std::endl;
//                return img::EasyImage();
                int nrIterations = configuration[figurename]["nrIterations"].as_int_or_default(1);
                figure = createMengerSponge(nrIterations);
            }

            Figures3D figures;
            if (fractal) {
                int nrIterations = configuration[figurename]["nrIterations"].as_int_or_default(1);
                double fractalScale = configuration[figurename]["fractalScale"].as_double_or_default(3.0);
                generateFractal(figure, figures, nrIterations, fractalScale);
            } else if (thickFigure) {
                double radius = configuration[figurename]["radius"].as_double_or_default(0.1);
                int n = configuration[figurename]["n"].as_int_or_default(180);
                int m = configuration[figurename]["m"].as_int_or_default(5);
                figure.generateThickFigure(figures, radius, n, m);
            }

            double rotateX = (configuration[figurename]["rotateX"].as_double_or_default(0)) * M_PI / 180;
            double rotateY = (configuration[figurename]["rotateY"].as_double_or_default(0)) * M_PI / 180;
            double rotateZ = (configuration[figurename]["rotateZ"].as_double_or_default(0)) * M_PI / 180;
            double scale = configuration[figurename]["scale"].as_double_or_default(1.0);
            std::vector<double> center = configuration[figurename]["center"].as_double_tuple_or_default({0, 0, 0});

            Matrix transformations = TransMatrix::scaleFigure(scale) *
                                     TransMatrix::rotateX(rotateX) *
                                     TransMatrix::rotateY(rotateY) *
                                     TransMatrix::rotateZ(rotateZ) *
                                     TransMatrix::translate(Vector3D::point(center[0], center[1], center[2])); // *
            // TransMatrix::eyePointTrans(Vector3D::point(eye_vector[0], eye_vector[1],
            //                                          eye_vector[2]));

            if (!fractal and !thickFigure) {
                figures.emplace_back(figure);
            }

            for (Figure figure1: figures) {
                if (!light)
                    figure1.ambientReflection = configuration[figurename]["color"].as_double_tuple_or_default(
                            {0, 0, 0});
                else
                    figure1.ambientReflection = configuration[figurename]["ambientReflection"].as_double_tuple_or_default(
                            {0, 0, 0});
                figure1.diffuseReflection = configuration[figurename]["diffuseReflection"].as_double_tuple_or_default(
                        {0, 0, 0});
                figure1.specularReflection = configuration[figurename]["specularReflection"].as_double_tuple_or_default(
                        {0, 0, 0});
                figure1.reflectionCoeff = configuration[figurename]["reflectionCoefficient"].as_double_or_default(0);
                figure1.applyTransformation(transformations);
                figures3D.emplace_back(figure1);
            }

        }

        Vector3D eyepoint = Vector3D::point(eye_vector[0], eye_vector[1], eye_vector[2]);
        Matrix eyepointM = TransMatrix::eyePointTrans(eyepoint);

        Lights3D lights;

        if (!light) {
            lights.emplace_back(new Light());
        } else {
            int nrLights = configuration["General"]["nrLights"].as_int_or_default(0);
            for (int i = 0; i < nrLights; ++i) {
                std::string figurename = "Light" + std::to_string(i);
                std::vector<double> ambient = configuration[figurename]["ambientLight"].as_double_tuple_or_default(
                        {0, 0, 0});
                std::vector<double> diffuse = configuration[figurename]["diffuseLight"].as_double_tuple_or_default(
                        {0, 0, 0});
                std::vector<double> specular = configuration[figurename]["specularLight"].as_double_tuple_or_default(
                        {0, 0, 0});

                if (configuration[figurename]["infinity"].as_bool_or_default(false)) {
                    std::vector<double> direction = configuration[figurename]["direction"].as_double_tuple_or_die();
                    Vector3D Ld = Vector3D::vector(direction[0], direction[1],
                                                   direction[2]); // + Vector3D::vector(.00, .00, 2); // kleine getallen zijn nodig om een vreemde bug op te lossen (die problemen gaf in de zbuffer als je er recht op keek)
                    Ld *= TransMatrix::eyePointTrans(Vector3D::point(eye_vector[0], eye_vector[1],
                                                                     eye_vector[2]));
                    Vector3D l = -Vector3D::normalise(Ld);
                    lights.emplace_back(new InfLight(ambient, diffuse, specular, l));
                } else if (configuration[figurename]["location"].exists()) {

                    std::vector<double> direction = configuration[figurename]["location"].as_double_tuple_or_die();
                    Vector3D pos = Vector3D::point(direction[0], direction[1], direction[2]);
                    int maskSize = configuration["General"]["shadowMask"].as_int_or_default(size);
                    if (!shadowEnabeled)
                        lights.emplace_back(new PointLight(ambient, diffuse, specular, pos, eyepoint));
                    else
                        lights.emplace_back(new PointLight(ambient, diffuse, specular, pos, maskSize, eyepoint));
                } else {
                    lights.emplace_back(new Light(ambient, diffuse, specular));
                }
            }
        }

        img::EasyImage image(size, size,
                             img::Color(backgroundcolor[0] * 255, backgroundcolor[1] * 255, backgroundcolor[2] * 255));
        bool background = configuration["General"]["background"].exists();
        std::string backgroundLocation;
        if (background) {
            backgroundLocation = configuration["General"]["background"].as_string_or_die();
        }

        Vector3D eyep = Vector3D::point(eye_vector[0], eye_vector[1], eye_vector[2]);
        if (type == "ZBufferedWireframe") {
            image.drawFigures3D(figures3D, size, true, eyep);
        } else if (type == "ZBuffering" or light) {
//            image.drawFigures3D(figures3D, size, true);
            image.drawZbufTriangFig(figures3D, size, lights,
                                    eyep, shadowEnabeled, background, backgroundLocation, false,
                                    std::__cxx11::string(), zBufferOutput);
            // TODO: fix draw function for zbuffered
        } else {
            image.drawFigures3D(figures3D, size, false, eyep);
        }
        return image;

    }


    return img::EasyImage();
}

int main(int argc, char const *argv[]) {

//    std::streambuf *psbuf, *backup;
//    std::ofstream filestr;
//    filestr.open("output.txt");
//
//    backup = std::cout.rdbuf();     // back up cout's streambuf
//
//    psbuf = filestr.rdbuf();        // get file's streambuf
//    std::cout.rdbuf(psbuf);         // assign streambuf to cout

    int retVal = 0;
    try {
        for (int i = 1; i < argc; ++i) {
            ini::Configuration conf;
            try {
                std::ifstream fin(argv[i]);
                fin >> conf;
                fin.close();
            }
            catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            img::EasyImage image = generate_image(conf);
            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string fileName(argv[i]);
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << argv[i] << std::endl;
            }
        }
    }
    catch (const std::bad_alloc &exception) {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }

//    std::cout.rdbuf(backup);        // restore cout's original streambuf
//
//    filestr.close();

    return retVal;
}
