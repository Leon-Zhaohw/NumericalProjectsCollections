/*
 *	email.h
 *
 *	Created by Ryoichi Ando on 8/4/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

#ifndef _EMAIL_H
#define _EMAIL_H

class email {
public:
	static void activate(); // Default deactivated, call here to activate
	static void print( const char *format, ...);
	static void setTitle( const char *title );
	static void send();
};

#endif