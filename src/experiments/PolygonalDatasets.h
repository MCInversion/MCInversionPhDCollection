
#include "pmp/Types.h"
#include <string>

/// \brief First polygonal example with disconnected bifurcation medial set
const std::vector unevenCrossPolyPts{
		pmp::Point2{39.507142, 14.544772},
		pmp::Point2{46.104261, 5.542495},
		pmp::Point2{61.36143, 4.906308},
		pmp::Point2{68.282948, 13.11281},
		pmp::Point2{66.153916, 31.426095},
		pmp::Point2{69.924933, 39.754365},
		pmp::Point2{111.082270, 39.723965},
		pmp::Point2{117.795930, 46.528945},
		pmp::Point2{117.765430, 66.419005},
		pmp::Point2{113.358230, 72.215125},
		pmp::Point2{89.030514, 71.743755},
		pmp::Point2{82.788235, 77.337235},
		pmp::Point2{87.565954, 122.613040},
		pmp::Point2{80.332640, 129.222760},
		pmp::Point2{68.372952, 128.366110},
		pmp::Point2{61.451434, 122.392570},
		pmp::Point2{57.089489, 85.394595},
		pmp::Point2{36.929786, 84.297475},
		pmp::Point2{10.835265, 83.544745},
		pmp::Point2{3.558908, 76.519305},
		pmp::Point2{3.558908, 57.450225},
		pmp::Point2{10.584357, 47.664785},
		pmp::Point2{32.413427, 48.166595},
		pmp::Point2{40.865615, 41.985195}
};

/// \brief Outer F outline for pair 0
const std::string svgPathPair00 = R"svg(
    <path
	   style="opacity:1;fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 41.99101,14.189914 6.597119,-9.0022801 15.257169,-0.63618 6.921518,8.2065001 -2.129032,18.31328 3.771017,8.32827 41.157339,-0.0304 6.71366,6.80498 -0.0305,19.89006 -4.4072,5.79612 -24.327718,-0.47137 -6.242279,5.59348 4.777719,45.275806 -7.233314,6.60972 -11.959688,-0.85665 -6.921518,-5.97354 -4.361945,-36.997966 -20.159699,-1.09713 -26.09452,-0.75272 -7.2763603,-7.02544 v -19.06908 l 7.0254503,-9.78545 21.82907,0.50182 8.452185,-6.1814 z"
       id="path288"
       sodipodi:nodetypes="ccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 0
const std::string svgPathPair01 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;opacity:1;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="M 48.967741,47.193549 40.451614,55 l 5.32258,8.516128 -0.354841,6.032259 -4.258063,9.225807 4.612904,2.838709 15.258063,0.354838 2.483871,-3.548388 2.483871,-7.806452 5.677419,-0.709676 8.16129,4.967743 5.32258,-2.129033 6.032259,-7.096774 -2.483871,-9.225806 -8.870968,-5.677421 -8.16129,1.774195 -9.225806,-2.483871 -2.129031,-9.225806 -7.096776,-6.03226 -4.612902,4.258067 z"
       id="path491" />
)svg";

/// \brief Outer F outline for pair 1
const std::string svgPathPair10 = R"svg(
    <path
	   style="fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 39.507141,13.835076 6.59712,-9.0022802 15.25717,-0.63618 6.92152,8.2065002 -2.12904,18.31328 3.77102,8.32827 41.157339,-0.0304 6.71366,6.80498 -0.0305,19.89006 -4.4072,5.79612 -24.327719,-0.47137 -6.24228,5.59348 4.77772,45.275804 -7.23331,6.60972 -11.95969,-0.85665 -6.92152,-5.97354 -4.36194,-36.997964 -20.1597,-1.09713 -26.09452,-0.75272 -7.27636,-7.02544 v -19.06908 l 7.02545,-9.78545 21.82907,0.50182 8.45218,-6.1814 z"
       id="path288"
       sodipodi:nodetypes="ccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 1
const std::string svgPathPair11 = R"svg(
    <path
    sstyle="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="M 93.217836,49.368896 81.85048,46.136053 71.576474,47.06779 65.946995,44.871651 61.409341,36.026613 53.38425,35.903399 l -5.21606,2.253249 -2.268553,4.343405 4.788501,7.850466 -1.345078,6.860967 -9.587311,2.255718 -3.133687,4.010275 0.227798,10.967142 3.049724,3.483688 18.216201,-1.288451 4.183297,1.094704 1.469706,10.858499 3.833505,4.230893 7.996524,-2.317974 0.917282,-6.474075 -1.124976,-12.214129 3.874897,-4.729449 9.432958,0.817044 7.92535,-4.892977 -1.919996,-6.055262 z"
       id="path491"
       sodipodi:nodetypes="cccccccccccccccccccccccccc" />
)svg";

/// \brief Outer F outline for pair 2
const std::string svgPathPair20 = R"svg(
    <path
       style="fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 34.238055,8.3150813 11.866206,-3.4822855 15.25717,-0.63618 15.70333,1.6828701 14.681853,3.0078401 11.298286,4.062819 12.80464,24.809577 1.94639,8.059524 -0.0305,19.89006 -2.90175,7.301573 -4.50592,9.815893 -6.49319,11.615292 -16.298619,27.461276 -7.23331,6.60972 -11.95969,-0.85665 -8.928791,-4.719 L 38.522236,111.78306 23.380713,102.65685 8.5770913,92.87141 3.558911,75.809616 V 56.740536 L 5.3152751,43.94418 16.104356,22.366021 l 7.950362,-8.43958 z"
       id="path288"
       sodipodi:nodetypes="ccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 2
const std::string svgPathPair21 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 90.379126,48.659219 -8.351227,-1.724779 -7.701425,-2.173102 -8.11335,-2.728397 -4.803783,-6.006328 -8.025091,-0.123214 -5.21606,2.253249 -2.268553,4.343405 0.885275,7.140789 -0.990239,6.328709 -6.038924,3.497653 -3.133687,4.010275 0.227798,10.967142 3.049724,3.483688 11.119427,-0.756193 8.530071,3.933414 4.219706,7.487531 3.833505,4.230893 7.996524,-2.317974 0.917282,-6.474075 -0.681428,-8.310903 5.205543,-6.681062 7.658764,-1.134569 7.92535,-4.892977 0.297746,-9.33752 z"
       id="path491"
       sodipodi:nodetypes="cccccccccccccccccccccccccc" />
)svg";

/// \brief Outer F outline for pair 3
const std::string svgPathPair30 = R"svg(
    <path
       style="fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 34.238055,8.3150813 11.866206,-3.4822855 15.25717,-0.63618 15.70333,1.6828701 14.681853,3.0078401 11.298286,4.062819 12.80464,24.809577 1.94639,8.059524 -0.0305,19.89006 -2.90175,7.301573 -4.50592,9.815893 -6.49319,11.615292 -16.298619,27.461276 -7.23331,6.60972 -11.95969,-0.85665 -8.928791,-4.719 L 38.522236,111.78306 23.380713,102.65685 8.5770913,92.87141 3.558911,75.809616 V 56.740536 L 5.3152751,43.94418 16.104356,22.366021 l 7.950362,-8.43958 z"
       id="path288"
       sodipodi:nodetypes="ccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 3
const std::string svgPathPair31 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 85.693545,36.282256 -4.080645,-1.596773 -18.80645,1.596773 -29.185484,14.10484 -1.508064,13.129032 4.524192,4.96774 6.653226,2.75 9.048388,-2.040321 10.82258,-4.080645 28.653226,-14.903225 -0.709679,-9.669357 z"
       id="path439" />
)svg";

/// \brief Outer F outline for pair 4
const std::string svgPathPair40 = R"svg(
    <path
       style="fill:none;stroke:#000000;stroke-width:1;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 39.507141,13.835076 6.59712,-9.0022802 15.25717,-0.63618 6.92152,8.2065002 -2.12904,18.31328 3.77102,8.32827 41.157339,-0.0304 6.71366,6.80498 -0.0305,19.89006 -4.4072,5.79612 -24.327719,-0.47137 -6.24228,5.59348 4.77772,45.275804 -7.23331,6.60972 -11.95969,-0.85665 -6.92152,-5.97354 -4.36194,-36.997964 -20.1597,-1.09713 -26.09452,-0.75272 -7.27636,-7.02544 v -19.06908 l 7.02545,-9.78545 21.82907,0.50182 8.45218,-6.1814 z"
       id="path288"
       sodipodi:nodetypes="ccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 4
const std::string svgPathPair41 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2.00000349,1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 84.995996,53.291363 13.45865,5.962993 3.593054,6.039291 8.22672,-1.561541 1.62942,-8.54683 -4.0938,-6.006278 -8.829882,-5.514636 -8.052867,0.39578 -6.002858,3.838865 0.07156,5.392356 z"
       id="path491"
       sodipodi:nodetypes="ccccccccccc" />
)svg";

/// \brief Inner G outline for pair 6
const std::string svgPathPair61 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 84.995996,53.291363 5.369238,7.52188 7.348753,4.392587 7.219163,0.965998 6.15672,-3.694266 0.81397,-7.292286 -4.0938,-6.006278 -8.829882,-5.514636 -8.052867,0.39578 -4.999223,3.211593 -1.06514,2.959144 z"
       id="path491"
       sodipodi:nodetypes="cccccccccccc" />
)svg";

/// \brief Inner G outline for pair 7
const std::string svgPathPair71 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 51.625119,68.722258 5.931384,8.346626 11.998499,1.272024 8.853994,-7.081536 0.374876,-10.052283 -4.470163,-9.393548 -7.073521,-6.392817 -8.244466,-1.137274 -6.206576,0.403964 -5.023675,3.085185 -1.606456,5.395494 0.949742,8.277808 z"
       id="path491"
       sodipodi:nodetypes="ccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 8
const std::string svgPathPair81 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="M 50.531056,72.928888 77.907178,63.983016 56.069915,42.629921 Z"
       id="path491"
       sodipodi:nodetypes="cccc" />
)svg";

/// \brief Inner G outline for pair 11
const std::string svgPathPair111 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 47.160089,68.848243 3.862953,0.53737 6.904906,0.6055 3.97766,1.891015 4.01837,0.02693 2.028355,-2.148588 4.99438,-7.466961 0.06656,-2.97427 -2.183816,-5.540785 -0.728726,-7.093821 -2.386188,-5.101784 -1.603755,-3.450252 -3.068597,-5.699643 0.651852,-4.079884 -0.56539,-6.24556 0.781439,-6.095647 -3.11224,-4.841568 L 58.0742,9.6474933 54.42049,9.5166176 51.811851,9.0089541 l -4.254087,3.7464289 -2.981742,4.452915 -1.298338,2.834554 1.019366,2.273668 1.683614,2.528925 -0.22957,2.769324 0.542811,2.811432 -1.441302,4.088867 0.368899,2.222201 0.803958,5.487846 -1.522021,4.528069 -0.676617,5.403279 -0.355749,4.703925 -0.195759,4.518771 0.261092,4.280321 z"
       id="path491"
       sodipodi:nodetypes="cccccccccccccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 12
const std::string svgPathPair121 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:1;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:2, 1;stroke-dashoffset:0;stroke-opacity:1"
       d="m 80.813598,55.416316 1.357848,2.334712 6.096176,5.561932 2.595817,1.186171 5.141077,0.982756 6.550784,0.06247 3.20412,-0.744193 2.53367,-1.821704 2.25734,-4.175562 0.4523,-4.775765 -0.31788,-4.264955 -3.52204,-5.115542 -2.54484,-1.550231 -7.221545,-0.3746 -5.362958,0.952074 -4.489806,1.791625 -2.831582,2.813897 -2.466187,3.534162 z"
       id="path491"
       sodipodi:nodetypes="ccccccccccccccccccc" />
)svg";

/// \brief Outer F outline for pair 16
const std::string svgPathPair160 = R"svg(
    <path
       style="fill:none;stroke:#000000;stroke-width:0.935438;stroke-linecap:square;stroke-linejoin:round;stroke-opacity:1"
       d="m 58.523369,23.447792 4.59688,-5.431163 11.047009,-4.182941 6.711603,6.046062 0.375758,4.106362 9.00908,8.449474 15.354491,-0.5716 11.09387,4.644115 5.71628,14.192127 -1.13073,17.833121 -4.92614,9.076575 -16.689457,0.476126 -5.582645,5.694954 -0.890733,17.412056 -5.526499,7.57003 -9.137616,-0.98112 -5.28828,-6.84141 L 67.00915,89.656269 57.231955,87.89175 46.919211,87.436062 41.35982,79.389917 v -21.83957 l 3.859623,-8.191019 8.561857,-2.847793 6.457758,-7.079477 z"
       id="path288"
       sodipodi:nodetypes="cccccccccccccccccccccccccc" />
)svg";

/// \brief Inner G outline for pair 16
const std::string svgPathPair161 = R"svg(
    <path
       style="fill:none;stroke:#601fb8;stroke-width:0.940893;stroke-linecap:round;stroke-linejoin:round;stroke-dasharray:1.88179, 0.940893;stroke-dashoffset:0;stroke-opacity:1"
       d="m 80.91692,54.334316 7.754203,0.218476 3.912064,6.02259 -1.853962,3.634412 1.733323,4.150095 6.797338,0.857937 4.628854,-4.712442 1.48074,-8.326073 -3.72025,-5.851141 -5.122035,-5.026524 -7.962978,-0.737882 -5.938798,1.579251 -3.159573,6.030842 z"
       id="path491"
       sodipodi:nodetypes="cccccccccccccc" />
)svg";

/// \brief An svg path parsing utility
[[nodiscard]] std::vector<pmp::Point2> ParsePolygonalSVGPath(const std::string& svgPath);