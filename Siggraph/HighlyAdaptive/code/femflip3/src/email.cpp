/*
 *	email.cpp
 *
 *	Created by Ryoichi Ando on 8/4/12
 *	Email: and@verygood.aid.design.kyushu-u.ac.jp
 *
 */

// See http://curl.haxx.se/libcurl/c/smtp-tls.html for detail

#include "util3.h"
#include "email.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <string>

#define ENABLE_LCURL		0
#if ENABLE_LCURL
#include <curl/curl.h>
#endif

static bool activated = false;
static int sent_count = 0;

/* This is a simple example showing how to send mail using libcurl's SMTP
 * capabilities. It builds on the simplesmtp.c example, adding some
 * authentication and transport security.
 */

#define FROM    "<youremail_address@server.com>"
#define TO      "<youremail_address@server.com>"

static std::string subject;
static std::string content;
static const char *payload_text[]={
	"To: " TO "\n",
	"From: " FROM "Program Notifier\n",
	NULL, //"Subject: SMTP TLS example message\n", // Slot: 2
	"\n", /* empty line to divide headers from body, see RFC5322 */
	NULL, // Slot: 4
	NULL
};

struct upload_status {
	int lines_read;
};

static size_t payload_source(void *ptr, size_t size, size_t nmemb, void *userp) {
	struct upload_status *upload_ctx = (struct upload_status *)userp;
	const char *data;
	
	if ((size == 0) || (nmemb == 0) || ((size*nmemb) < 1)) {
		return 0;
	}
	
	data = payload_text[upload_ctx->lines_read];
	
	if (data) {
		size_t len = strlen(data);
		memcpy(ptr, data, len);
		upload_ctx->lines_read ++;
		return len;
	}
	return 0;
}

void email::activate() {
	activated = true;
}

void email::print( const char *format, ...) {
	// Append content
	char message[512];
	va_list args;
	va_start(args, format);
	vsprintf(message, format, args);
	va_end(args);
	content += message;
}

void email::setTitle( const char *title ) {
	subject = "Subject: ";
	subject += title;
}

void email::send() {
	if( ! activated ) {
		printf( "Send mail attempted: %s", content.c_str() );
		exit(0);
	} else {
		// Check maximum sent count
		if( sent_count > 3 ) {
			dump( "Too much email sent. Exiting for safety...\n" );
			exit(0);
		} else {
			dump( "Sending an email...\n" );
			dump( "%s\n" , subject.c_str());
			dump( "Content: %s\n", content.c_str() );
		}
		
		// Correct return code
		for( uint n=0; n<content.size(); n++ ) {
			if( content[n] == '\n' && content[n+1] == 0 ) content[n] = 0;
		}
		// Set message
		payload_text[2] = subject.c_str();
		payload_text[4] = content.c_str();
	}
#if ENABLE_LCURL
	CURL *curl;
	CURLcode res;
	struct curl_slist *recipients = NULL;
	struct upload_status upload_ctx;
	upload_ctx.lines_read = 0;
	
	curl = curl_easy_init();
	if (curl) {
		/* This is the URL for your mailserver. Note the use of port 587 here,
		 * instead of the normal SMTP port (25). Port 587 is commonly used for
		 * secure mail submission (see RFC4403), but you should use whatever
		 * matches your server configuration. */
		curl_easy_setopt(curl, CURLOPT_URL, "smtp://smtp.yoursever.com:587");
		
		/* In this example, we'll start with a plain text connection, and upgrade
		 * to Transport Layer Security (TLS) using the STARTTLS command. Be careful
		 * of using CURLUSESSL_TRY here, because if TLS upgrade fails, the transfer
		 * will continue anyway - see the security discussion in the libcurl
		 * tutorial for more details. */
		curl_easy_setopt(curl, CURLOPT_USE_SSL, (long)CURLUSESSL_ALL);
		
		/* If your server doesn't have a valid certificate, then you can disable
		 * part of the Transport Layer Security protection by setting the
		 * CURLOPT_SSL_VERIFYPEER and CURLOPT_SSL_VERIFYHOST options to 0 (false).
		 *   curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
		 *   curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);
		 * That is, in general, a bad idea. It is still better than sending your
		 * authentication details in plain text though.
		 * Instead, you should get the issuer certificate (or the host certificate
		 * if the certificate is self-signed) and add it to the set of certificates
		 * that are known to libcurl using CURLOPT_CAINFO and/or CURLOPT_CAPATH. See
		 * docs/SSLCERTS for more information.
		 */
		// curl_easy_setopt(curl, CURLOPT_CAINFO, "/path/to/certificate.pem");
		
		/* A common reason for requiring transport security is to protect
		 * authentication details (user names and passwords) from being "snooped"
		 * on the network. Here is how the user name and password are provided: */
		curl_easy_setopt(curl, CURLOPT_USERNAME, "email_address@yourserver.com");
		curl_easy_setopt(curl, CURLOPT_PASSWORD, "yourpassword");
		
		/* value for envelope reverse-path */
		curl_easy_setopt(curl, CURLOPT_MAIL_FROM, FROM);
		/* Add two recipients, in this particular case they correspond to the
		 * To: and Cc: addressees in the header, but they could be any kind of
		 * recipient. */
		recipients = curl_slist_append(recipients, TO);
		curl_easy_setopt(curl, CURLOPT_MAIL_RCPT, recipients);
		
		/* In this case, we're using a callback function to specify the data. You
		 * could just use the CURLOPT_READDATA option to specify a FILE pointer to
		 * read from.
		 */
		curl_easy_setopt(curl, CURLOPT_READFUNCTION, payload_source);
		curl_easy_setopt(curl, CURLOPT_READDATA, &upload_ctx);
		
		/* Since the traffic will be encrypted, it is very useful to turn on debug
		 * information within libcurl to see what is happening during the transfer.
		 */
		// curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
		curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
		
		/* send the message (including headers) */
		res = curl_easy_perform(curl);
		/* Check for errors */
		if(res != CURLE_OK) {
			fprintf(stderr, "curl_easy_perform() failed: %s\n",
					curl_easy_strerror(res));
			exit(0);
		}
		/* free the list of recipients and clean up */
		curl_slist_free_all(recipients);
		curl_easy_cleanup(curl);
		
		// Clear message body
		content.clear();
		
		// Record sent count
		sent_count ++;
	}
#endif
}
