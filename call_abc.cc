#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

// procedures to start and stop the ABC framework
// (should be called before and after the ABC procedures are called)
extern "C"
{
extern void   Abc_Start();
extern void   Abc_Stop();

// procedures to get the ABC framework and execute commands in it
extern void * Abc_FrameGetGlobalFrame();
extern int    Cmd_CommandExecute( void * pAbc, char * sCommand );
}

int call_abc(int flag)
{

    // variables
    void * pAbc;
    FILE *fp;
    char Command1[500], Command2[500], Command3[500], Command4[500], Command5[500], Command6[500], Command7[500];
    
     //////////////////////////////////////////////////////////////////////////
     // start the ABC framework
     Abc_Start();
     pAbc = Abc_FrameGetGlobalFrame();
     
    if(flag == 0)  
    {
	    sprintf( Command1, "read_blif ./blif_files/ckt_org_simu.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "write_verilog ./verilog_files/ckt_org_simu.v");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }
    else if(flag == 1)   
    {
        sprintf( Command1, "read_blif ./blif_files/ckt_org_sim_simu.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "sweep; write_verilog ./verilog_files/ckt_org_sim_simu.v");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }
    else if(flag == 2)   
    {
        sprintf( Command1, "read_blif ./blif_files/ckt_assure_simu.blif");
		if ( Cmd_CommandExecute( pAbc, Command1 ) )
		{
	    	fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	    	return 1;
		} 
	    
		sprintf( Command1, "write_verilog ./verilog_files/ckt_assure_simu.v");
		if ( Cmd_CommandExecute( pAbc, Command1 ) )
		{
	    	fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	    	return 1;
		}
    }
    else if(flag == 3)
    {
	    sprintf( Command1, "read_blif ./blif_files/ckt_org_simu_copy.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "write_verilog ./verilog_files/ckt_org_simu_copy.v");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }
    else if(flag == 4)   
    {
        sprintf( Command1, "read_blif ./blif_files/ckt_org_sim_simu_copy.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "sweep; write_verilog ./verilog_files/ckt_org_sim_simu_copy.v");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }
    else if(flag == 5)   
    {
        sprintf( Command1, "read_blif ./blif_files/ckt_org_sim_simu_copy.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "sweep; write_verilog ./verilog_files/ckt_org_sim_simu_copy.v");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }
    else if(flag == 6)   
    {
        sprintf( Command1, "read_blif ./blif_files/compare_ckt_max_sat.blif");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
	    
	    sprintf( Command1, "sweep;strash; write_cnf compare_ckt_max_sat.cnf");
	    if ( Cmd_CommandExecute( pAbc, Command1 ) )
	    {
	        fprintf( stdout, "Cannot execute command \"%s\".\n", Command1 );
	        return 1;
	    }
    }

    //////////////////////////////////////////////////////////////////////////
    // stop the ABC framework
    Abc_Stop();
    return 0;   
}

