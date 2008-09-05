#include <gtest/gtest.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include "ParmParser.h"

using namespace std;

namespace {
    
// Utility function to construct an argument list (argc, argv) from
// a given string
void GetArgList(string s, int& argc, char**& argv ) {
    // First, count number of arguments
    istringstream in(s);
    argc = 0;
    string str;
    while ( in >> str ) {
        ++argc;
    }
    
    // Allocate memory for argv
    argv = (char **) malloc( argc * sizeof(char *) );

    // Reset input stream and get arguments
    in.clear();
    in.str(s);
    for (int i=0; i<argc; ++i) {
        in >> str;
        argv[i] = (char *) malloc ( (str.length()+1) * sizeof(char) );
        strcpy( argv[i], str.c_str() );
    }
}

ParmParser* GetParser(string s) {
    int argc;
    char **argv;
    GetArgList( s, argc, argv );
    ParmParser* parser = new ParmParser( argc, argv );
    return parser;
}

TEST( ParmParserTest, True ) {
    EXPECT_EQ(1, 1);
}

TEST( ParmParserTest, CorrectArguments ) {
    string cmd = "a.out -b 1 -flag -fname myFile -ny 20 -nx 10 -length 3.14";
    ParmParser* parser = GetParser( cmd );
    int nx = parser->getInt( "nx", "description of nx", 100 );
    EXPECT_EQ( 10, nx );

    int ny = parser->getInt( "ny", "description of ny", 200 );
    EXPECT_EQ( 20, ny );

    double length = parser->getDouble( "length", "description", 10. );
    EXPECT_DOUBLE_EQ( 3.14, length );

    string fname = parser->getString( "fname", "description", "name" );
    EXPECT_EQ( "myFile", fname );
    
    bool on = parser->getBool( "b", "description", 0 );
    EXPECT_EQ( true, on );

    bool flag = parser->getFlag( "flag", "description" );
    EXPECT_EQ( true, flag );
    
    bool flag2 = parser->getFlag( "flag2", "description" );
    EXPECT_EQ( false, flag2 );

    EXPECT_EQ( true, parser->inputIsValid() );
    EXPECT_EQ( "a.out -nx 10 -ny 20 -length 3.14 -fname myFile -b 1 -flag",
        parser->getParameters() );
    
    delete parser;
}

TEST( ParmParserTest, BadArg ) {
    string cmd = "a.out -nx 10 bad";
    ParmParser* parser = GetParser( cmd );
    int nx = parser->getInt( "nx", "description", 100 );
    EXPECT_EQ( 10, nx );
    EXPECT_EQ( false, parser->inputIsValid() );
    EXPECT_EQ( "a.out -nx 10", parser->getParameters() );
    delete parser;
}

TEST( ParmParserTest, MissingValueAtEnd ) {
    string cmd = "a.out -nx";
    ParmParser* parser = GetParser( cmd );
    int nx = parser->getInt( "nx", "desc", 100 );
    EXPECT_EQ( 100, nx );
    EXPECT_EQ( false, parser->inputIsValid() );
    EXPECT_EQ( "a.out", parser->getParameters() );
    delete parser;
}

TEST( ParmParserTest, MissingValueInMiddle ) {
    string cmd = "a.out -nx -ny 20";
    ParmParser* parser = GetParser( cmd );
    int nx = parser->getInt( "nx", "desc", 100 );
    EXPECT_EQ( 100, nx );
    int ny = parser->getInt( "ny", "desc", 200 );
    EXPECT_EQ( 20, ny );
    EXPECT_EQ( false, parser->inputIsValid() );
    EXPECT_EQ( "a.out -ny 20", parser->getParameters() );
    delete parser;
}

} // namelist