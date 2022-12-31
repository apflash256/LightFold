#include <iostream>
#include <imageio.h>
#include <math/vecmathalt.h>

int main(void) {
	std::cout << "Hello World!";
	float testimg[30000];
	for (int i = 0; i < 30000; i++) {
		testimg[i] = 1.0f * i;
	}
	char testfname[] = "testimg.exr";
	int a = WriteEXR(testimg, 100, 100, testfname);
	return a;
}