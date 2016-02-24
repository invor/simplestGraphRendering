#version 130

flat in int id;

out int fragColour;

void main()
{
    fragColour = id;
}