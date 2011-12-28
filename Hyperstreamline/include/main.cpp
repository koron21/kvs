#include "CellLocatorBIH.h"
#include <kvs/Timer>
#include <kvs/UnstructuredVolumeImporter>

//#define WIN_4131
#define WIN_HOME

using namespace std;

int main(int argc, char** argv)
{
	string filename = (argc == 2) ? argv[1] : 
#if defined WIN_HOME
		"D:\\Koron\\Work\\Data\\KVSML\\UnstructuredVolumeObject\\engine\\v6engine.kvsml";
#elif defined WIN_4131
        "C:\\Users\\Kaku\\Dropbox\\Work\\Viz\\Hyper-Streamline\\data\\engine\\v6engine.kvsml";
#elif defined LINUX
		"/home/koron/Dropbox/Work/Viz/Hyper-Streamline/data/engine/v6engine.kvsml";
#endif
	kvs::UnstructuredVolumeObject* volume = new kvs::UnstructuredVolumeImporter(filename);

   	//single threaded mode
	cout << "Single threaded mode\n";
	kvs::CellLocatorBIH* locator_s = new kvs::CellLocatorBIH(volume);
	locator_s->build();


	//double threaded mode
	cout << "Double threaded mode\n";
	kvs::CellLocatorBIH* locator_d = new kvs::CellLocatorBIH(volume);
	locator_d->setParallel();
	locator_d->build();
	
	delete locator_s;
	delete locator_d;
	return 0;
}