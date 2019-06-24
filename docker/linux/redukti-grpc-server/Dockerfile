### Runs the Redukti GRPC Server
FROM redukti-ubuntu-risk-base:latest

LABEL description="redukti-grpc-server" maintainer="redukti@gmail.com"

EXPOSE 9001

CMD /Software/redukti/bin/reduktiserv -s /Software/redukti/scripts/pricing.lua -p 9001
