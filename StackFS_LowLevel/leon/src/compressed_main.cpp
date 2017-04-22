/*****************************************************************************
 *   Leon: reference free compression for NGS reads
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Benoit, G.Rizk, C.Lemaitre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <gatb/gatb_core.hpp>
#include <Leon.hpp>


using namespace std;


/********************************************************************************/


void displayVersion(std::ostream& os){

	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	os << "* Leon version "<< LEON_VERSION_MAJOR << "."
		<< LEON_VERSION_MINOR << "."
		<< LEON_VERSION_PATCH
		<< "                      *" << endl; //<< " AGPL licence" <<endl;
	os << "* Using gatb-core version "<< STR_LIBRARY_VERSION <<  "           *" << endl;
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
}

int main (int argc, char* argv[])
{


	if(argc > 1 && (   strcmp(argv[1],STR_VERSION)==0 || strcmp(argv[1],"-v")==0    )     ){
		displayVersion(cout);
		return EXIT_FAILURE;
	}
	// We define a try/catch block in case some method fails
	try
	{
		Leon *leon = new Leon();	
		leon->run (argc, argv);
		string dir = System::file().getDirectory(leon->_inputFilename);
		int block_count = 0;
		bool isFasta = true;
		IBank* whole_bank = Bank::open(leon->_inputFilename);
		int64_t seqCount = whole_bank->estimateNbItems();
		string temp_file(".tmp_file_compress");
		if(leon->_inputFilename.find(".fq") !=  string::npos || leon->_inputFilename.find(".fastq") !=  string::npos){
			if(! leon->getParser()->saw (Leon::STR_DNA_ONLY) && ! leon->getParser()->saw (Leon::STR_NOQUAL)){
				isFasta = false;
			}
			temp_file += ".fq";
		}else
			temp_file += ".fa";
		temp_file = dir+"/"+temp_file;
		Iterator<Sequence>* itSeq = leon->createIterator<Sequence> (whole_bank->iterator(),seqCount,
				"Creating Temp Files"
				);
		ifstream fin;
		fin.open(leon->_inputFilename.c_str());
		string line;
		leon->orig_block_size->clear();
		leon->seq_per_block->clear();
		bool more_lines = true;
		bool first_time = true;
		while(more_lines){
			Leon *leon1  = new Leon();
			leon1->run (argc, argv);
			ofstream fout;
			int j = 0, size = 0;
			fout.open(temp_file.c_str());   //create a new file for each block
			if (!fout.good()){
				cerr << "I/O error while reading file." << endl;
			}
			if(!first_time){                //do not execute the first time
				j++;
				fout<<line.c_str()<<endl;
				size += line.size()+1;
			}
			while(getline(fin, line))
			{
				if(isFasta && line.at(0)=='>')
					j++;
				else if(!isFasta && line.at(0)=='@')
					j++;
				if(j>leon->READ_PER_BLOCK){
					j--;
					break;
				}
				fout<<line.c_str()<<endl;
				size += line.size()+1;  //+1 for the end of line character
			}
			first_time = false;
			if(j==0)
				break;
			if(j<leon->READ_PER_BLOCK){
				more_lines = false;
			}
			leon->orig_block_size->push_back(size);
			leon->seq_per_block->push_back(j);
			string file_name (leon->_base_outputFilename);
			file_name += "_"+to_string(block_count)+".leon";
			leon->outputFileNames->push_back(file_name);
			fout.close();
			leon1->executeCompression(block_count, temp_file.c_str());
			block_count++;
		}
		leon->saveConfig();
	}
	catch (gatb::core::system::Exception& e)
	{

		cerr << "EXCEPTION: " << e.getMessage() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

