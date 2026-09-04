# Multi-stage production container for Dark Calculator Suite
FROM node:22-alpine AS runtime

WORKDIR /app

# Copy dependency manifests
COPY package*.json ./

# Copy application source code
COPY . .

# Environment settings
ENV NODE_ENV=production
ENV PORT=3000

# Expose standard application port
EXPOSE 3000

# Healthcheck
HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
  CMD node -e "require('http').get('http://localhost:3000/api/health', (r) => process.exit(r.statusCode === 200 ? 0 : 1))"

# Start application server
CMD ["npm", "start"]
