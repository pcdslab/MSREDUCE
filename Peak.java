/*
This program accepts a folder of .dta files, percentage data 
to be retained and the output folder as arguments. The output 
folder contains .ms2 file which can be directly used with 
Crux-Tide search software.

java msreduce ./inputFolder xx ./outputFolder
xx = percentage data to be retained e.g. 10, 20, 30.

Copyright (C) Muaaz Gul Awan and Fahad Saeed  name of author

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/



public class Peak {
    
    public float intensity;
    public float m_z;
            
}
