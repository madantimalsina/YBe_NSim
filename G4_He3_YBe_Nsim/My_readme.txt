// This is for the Geometry based on our He3 tubes set up at SDSMT
 Folder G4_He3_SDSMT_setup:

E) Neutron detected for the SDSMT setup
New on the folder MyG4_SDSMT_setup_He3

6). Three Teflon layer: one on the the back and two on the side 
 (Calorimeter size make bigger to track all the information).
7). Changed all the name B4a into He3!

From old folder MyG4_He3fill_3Tubes_good_basic2
1). Source holder filled with Plexiglass.
2). Teflon as a Gap materials.
3). Outer layer of He3 tube made with Al (not final yet)
4). Inner layer of He3 tube filled with He3 gas
5). The Teflon or any layer below the He3 tubes.


// This is for the Geometry based on our He3 tubes set up
Same as random geometry below C.1 and C.2 but for the our Geometry 
D) Neutron detected from individual tubes (Including outer layer of 			   	Tubes in the our Geometry)
MyG4_He3fill_3Tubes_good_basic
1). Source holder filled with Plexiglass.
2). Teflon as a Gap materials.
3). Outer layer of He3 tube made with Al (not final yet)
4). Inner layer of He3 tube filled with He3 gas

See on MyG4_He3fill_3Tubes_good_basic2 for all including 
5). The Teflon or any layer behind the He3 tubes.

// This is for random geometry I tired
Now in the folder My_He3fill_2nTry
C) Exciting so far -->Neutron detected from individual tubes (Including outer layer of 			 Tubes in the Geometry)

Folder --> My_He3fill_3Tubes_good
1) So far Geometry is almost complete ( Source holder --> Teflon --> outer layer of He3 Tube --> Inner layer of He3 Tube fill with He3 gas).
2) Neutron recorded in the different tubes will save in the individual Ntuples. We can get number of neutron recorded in each tube as well as in the all tubes. 

// B) -->Neutron detected from individual tubes (But there is no outer layer of Tubes in 	 the Geometry)
After A.2
1) The folder My_He3fill_3Tubes_work:
I recorded the neutron numbers from individual tubes and save the data in a different Ntuples and histogram so that I will have number of neutron recorded for different tubes.
Once I have data for the different different tube I can get the total number by adding them!  


//////////
It is old all in folder My_He3fill_1stTry
A) -->Working on Geometry 
1). Here up to B4a_My_He3fill_work3 I am fixing the geometry. So far it works for the different thickness. But I am not sure all the three tube are taking data or not. So my next step is test if all the tubes are working and change based on the Andrew suggestion.

2).At B4a_My_He3fill_work4 I tried to check if the individual tube are tracking energy correctly or not. It's give me some wired output: all the tubes 2 sides and 1 center (when commented with /**/) give me the different number of neutron which is expected but sum is not equal to the neutron from total number of tubes. Neutron form all tubes is less then neutron from center tube only! Let's try to record ntuples from individual tubes.
  