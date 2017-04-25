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
		bool isFasta = true, noHeader = false;
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
		if(leon->getParser()->saw (Leon::STR_NOHEADER))
			noHeader = true;
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
		itSeq->first();
		while(! itSeq->isDone()){
			Leon *leon1  = new Leon();
			leon1->run (argc, argv);
			ofstream fout;
			int j = 0, size = 0, readid = 0;
			bool  reading = true;
			fout.open(temp_file.c_str());   //create a new file for each block
			if (!fout.good()){
				cerr << "I/O error while reading file." << endl;
			}
			while(! itSeq->isDone() && j < leon->READ_PER_BLOCK)
			{
                                stringstream sint;
                                sint << readid;
                                if(!noHeader)
				{
					string line = itSeq->item().getComment();
					if(isFasta)
						fout<<">";
					else
						fout<<"@";
					fout<<line<<'\n';
					size += 1+line.size()+1;
				}
				else
				{
                                        if(isFasta)
                                                fout<< "> "<<sint.str()<<'\n';
                                        else
                                                fout<< "@ "<<sint.str()<<'\n';
                                        readid++;
					size+= 2+sint.str().size()+1;
                                }
                                fout<< itSeq->item().getDataBuffer()<<'\n';
				size+= strlen(itSeq->item().getDataBuffer())+1;
				if( !isFasta)
                                {	
					string line = itSeq->item().getQuality();                                       
                                       	fout<<"+\n";
                                        fout<<line<<'\n';
					size+= 2+line.size()+1;                                        
                                }
				j++;
				itSeq->next();
			}
			leon->orig_block_size->push_back(size);
			leon->seq_per_block->push_back(j);
			fout.close();
			leon1->executeCompression(block_count, temp_file.c_str());
			block_count++;
		}
		//leon->saveConfig();
	}
	catch (gatb::core::system::Exception& e)
	{

		cerr << "EXCEPTION: " << e.getMessage() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

